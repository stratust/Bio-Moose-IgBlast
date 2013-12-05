use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::IgBlast {
    use MooseX::StrictConstructor;
    has 'molecule'              => ( is => 'ro', isa => 'Str',  required => 1 );
    has 'version'               => ( is => 'ro', isa => 'Str',  required => 1 );
    has 'query_id'              => ( is => 'ro', isa => 'Str',  required => 1 );
    has 'database'              => ( is => 'ro', isa => 'Str',  required => 1 );
    has 'domain_classification' => ( is => 'ro', isa => 'Str',  required => 0 );
    has 'converted_sequence'    => ( is => 'rw', isa => 'Bool', required => 1, default => 0 );
    has 'rearrangement_summary' => (
        is       => 'ro',
        isa      => 'Bio::Moose::IgBlast::Rearrangement',
        required => 0
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
                if ( $self->rearrangement_summary->top_V_match =~ /IGH/i ) {

                    if ( $query->sequence =~ /($r)(\S+)TGGGGC/i ) {
                        $cdr3_seq = $2 . "|";
                    }
                    elsif ( $query->sequence =~ /($r)(\S+)/i ) {
                        $cdr3_seq = $2;
                    }
                }

                # If light chain
                elsif ( $self->rearrangement_summary->top_V_match =~ /IG[LK]/i ) {

                    if ( $query->sequence =~ /($r)(\S+)TCCTGT/i ) {
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

                if ( $self->rearrangement_summary->top_V_match =~ /IGH/i ) {
                    if ( $query->translation =~ /($r)(\S+)WG/i ) {
                        $cdr3_seq = $2 . "|";
                    }
                    elsif ( $query->translation =~ /($r)(\S+)/i ) {
                        $cdr3_seq = $2;
                    }
                }

                # If light chain
                elsif ( $self->rearrangement_summary->top_V_match =~ /IG[LK]/i ) {
                    if ( $query->translation =~ /($r)(\S+)SC/i ) {
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

}

 # ABSTRACT: turns baubles into trinkets
