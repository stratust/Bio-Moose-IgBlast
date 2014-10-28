use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::IgBlast {
    use MooseX::StrictConstructor;
    use Text::Brew qw(distance);

    has 'molecule'              => ( is => 'ro', isa => 'Str',  required => 1 );
    has 'version'               => ( is => 'ro', isa => 'Str',  required => 1 );
    has 'query_id'              => ( is => 'ro', isa => 'Str',  required => 1 );
    has 'query_length'          => ( is => 'ro', isa => 'Int',  required => 0 );
    has 'database'              => ( is => 'ro', isa => 'Str',  required => 1 );
    has 'domain_classification' => ( is => 'ro', isa => 'Str',  required => 0 );
    has 'converted_sequence'    => ( is => 'rw', isa => 'Bool', required => 1, default => 0 );
    has 'rearrangement_summary' => (
        is       => 'ro',
        isa      => 'Bio::Moose::IgBlast::Rearrangement',
        required => 0
    );
    
    has 'germline_score' => (
        is            => 'ro',
        isa           => 'ArrayRef[Bio::Moose::IgBlast::GermlineScore]',
        required      => 0,
        traits        => ['Array'],
        documentation => 'ArrayRef of germline scores',
        handles       => {
            all_scores   => 'elements',
            add_score    => 'push',
            next_score   => 'shift',
            count_scores => 'count',
        },
    );

    has 'junction_details' => (
        is       => 'ro',
        isa      => 'Bio::Moose::IgBlast::Junction',
        required => 0
    );
    has 'alignments' => (
        is       => 'ro',
        isa      => 'Bio::Moose::IgBlast::Alignment',
        required => 0
    );
    has 'hit_table' => (
        is       => 'ro',
        isa      => 'Bio::Moose::IgBlast::HitTable',
        required => 0
    );

    has 'rendered_alignment' => (
        is       => 'ro',
        isa      => 'Bio::Moose::IgBlast::RenderedAlignment',
        required => 0
    );

    has 'chain_type' => (
        is      => 'ro',
        isa     => 'Str',
        lazy    => 1,
        builder => '_infer_chain_type',
    );
 
    has 'best_V' => (
        is      => 'ro',
        isa     => 'Str',
        lazy    => 1,
        builder => '_infer_best_V',
    );
   
    has 'best_D' => (
        is      => 'ro',
        isa     => 'Str',
        lazy    => 1,
        builder => '_infer_best_D',
    );

    has 'best_J' => (
        is      => 'ro',
        isa     => 'Str',
        lazy    => 1,
        builder => '_infer_best_J',
    );

    # store any information you want in this part
    has misc     => ( is => 'rw', isa => 'Any' );
    has init_pos => ( is => 'ro', isa => 'Int' );


    method infer_CDR3_nt {
        my $cdr3_seq = "N/A";
    
        if ( $self->rendered_alignment ) {
            my $query = $self->rendered_alignment->query;

            if ( $query->sub_regions_sequence->has_FWR3 ) {

                my $r = quotemeta $query->sub_regions_sequence->FWR3;

                # If heavy chain
                if ( $self->rearrangement_summary->chain_type =~ /VH/i ) {

                    if ( $query->sequence =~ /($r)(\S+)TGGGG[ATCG]/i ) {
                        $cdr3_seq = $2 . "|";
                    }
                    elsif ( $query->sequence =~ /($r)(\S+)/i ) {
                        $cdr3_seq = $2;
                    }
                }

                # If light chain
                elsif ( $self->rearrangement_summary->chain_type =~ /V[LK]/i ) {

                    if ( $query->sequence =~ /($r)(\S+)TT[CT]GG[TC]/i ) {
                        $cdr3_seq = $2 . "|";
                    }
                    elsif ( $query->sequence =~ /($r)(\S+)/i ) {
                        $cdr3_seq = $2;
                    }
                }
            }
        }
        return $cdr3_seq;
    }


    method infer_CDR3_nt_length {
        my $cdr3_seq = $self->infer_CDR3_nt;
        my $len = "N/A";
        if ($cdr3_seq ne "N/A"){
            $cdr3_seq =~ s/\|//g;
            $len = length $cdr3_seq;
        }
        return $len;
    }


    method infer_CDR3_aa {
        my $cdr3_seq = "N/A";

        if ( $self->rendered_alignment ) {
            my $query = $self->rendered_alignment->query;

            if ( $query->sub_regions_translation->has_FWR3 ) {

                my $r = quotemeta $query->sub_regions_translation->FWR3;

                if ( $self->rearrangement_summary->chain_type =~ /VH/i ) {
                    if ( $query->translation =~ /($r)(\S+)WG/i ) {
                        $cdr3_seq = $2 . "|";
                    }
                    elsif ( $query->translation =~ /($r)(\S+)/i ) {
                        $cdr3_seq = $2;
                    }
                }

                # If light chain
                elsif ( $self->rearrangement_summary->chain_type =~ /V[LK]/i ) {
                    if ( $query->translation =~ /($r)(\S+)FG/i ) {
                        $cdr3_seq = $2 . "|";
                    }
                    elsif ( $query->translation =~ /($r)(\S+)/i ) {
                        $cdr3_seq = $2;
                    }
                }
            }
        }
        return $cdr3_seq;
    }


    method infer_CDR3_aa_length {
        my $cdr3_seq = $self->infer_CDR3_aa;
        my $length = "N/A";
        if ($cdr3_seq ne "N/A"){
            $cdr3_seq =~ s/\|//g;
            $length = length $cdr3_seq;
        }
        return $length;
    }


    method infer_aa_diff (Str $region where [qr/FWR[123]/,qr/CDR[12]/]) {
        my ( $mismatches, $insertions, $deletions ) = ('N/A') x 3;

        if ( $self->rendered_alignment ) {

            if ( $self->rearrangement_summary->stop_codon =~ /no/i ) {

                my $query  = $self->rendered_alignment->query;
                my $best_V = $self->rendered_alignment->best_V;
                my $r      = "has_$region";

                if ( $query->sub_regions_translation->$r && $best_V->sub_regions_translation->$r ) {

                    my ( $distance, $arrayref_edits ) = distance(
                        $best_V->sub_regions_translation->$region,
                        $query->sub_regions_translation->$region
                    );
                    my @subs    = grep( /SUBST/, @{$arrayref_edits} );
                    my @inserts = grep( /DEL/,   @{$arrayref_edits} );
                    my @dels    = grep( /INS/,   @{$arrayref_edits} );

                    $mismatches = scalar @subs;
                    $insertions = scalar @inserts;
                    $deletions  = scalar @dels;
                }
            }
        }
        return ( $mismatches, $insertions, $deletions );
    }

    
    method _infer_chain_type {
        my $type = "N/A";
        if ( $self->rearrangement_summary ) {
            my $r = $self->rearrangement_summary;
            if ( $r->top_V_match ) {
                if ( $r->top_V_match =~ /IGH/i ) {
                    $type = 'gamma';
                }
                elsif ($r->top_V_match =~ /IGL/i){
                    $type = 'lambda';
                }
                elsif ($r->top_V_match =~ /IGK/i){
                    $type = 'kappa';
                } 

            }
        }
        return $type;
    }


    # Reliable means E value less than 1
    method is_reliable {
        my $answer = 0;
        if ($self->germline_score){
            my $germ = $self->germline_score->[0];
            $answer = 1 if $germ->evalue < 1;
        } 
        return $answer;
    }


    # complete means all V[D]J regions
    method is_complete {
        my $answer = 0;
        if ( $self->rearrangement_summary ) {
            my $r = $self->rearrangement_summary;
            if ( $r->top_V_match && $r->top_J_match ) {
                if ( $r->top_V_match ne 'N/A' && $r->top_J_match ne 'N/A' ) {
                    if ( $self->chain_type =~ /kappa|lambda/i ) {
                        $answer = 1;
                    }
                    elsif ( $self->chain_type =~ /gamma/i ) {
                        $answer = 1 if $r->top_D_match && $r->top_D_match ne 'N/A';
                    }
                }
            }
        }

        return $answer;
    }


    method _infer_best_V {
        my $chain = 'N/A';
        if ($self->is_reliable){
             my $r = $self->rearrangement_summary;
             $chain = $r->top_V_match;
        }
        return $chain;
    }


    method _infer_best_D {
        my $chain = 'N/A';
        if ($self->chain_type =~ /gamma/i){
             my $r = $self->rearrangement_summary;
             $chain = $r->top_D_match if $r->top_D_match;
        }
        return $chain;
    }


    method _infer_best_J {
        my $chain = 'N/A';
        if ( $self->rearrangement_summary ) {
            my $r = $self->rearrangement_summary;
            $chain = $r->top_J_match if $r->top_J_match;
        }
        return $chain;
    }

}

 # ABSTRACT: turns baubles into trinkets
