use MooseX::Declare;
use feature 'say';
use Method::Signatures::Modifiers;

class Bio::Moose::IgBlast {
    use MooseX::StrictConstructor;
    use Text::Brew qw(distance);
    use Bio::Seq;
    use Data::Printer;

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
    has igblast_block => ( is => 'ro', isa => 'Str' );

    has 'mismatches' => (
        is      => 'rw',
        isa     => 'HashRef|Undef',
        lazy    => 1,
        builder => '_build_mismatches',
    );

    method _build_mismatches {
        my $mismatches;

        my ( $germ_nt, $germ_aa, $query_nt, $query_aa );

        # Check germline best V region;
        if ( $self->rendered_alignment && $self->rendered_alignment->best_V ) {
            my $V = $self->rendered_alignment->best_V;
            if ( $V->sub_regions_sequence && $V->sub_regions_translation ) {
                $germ_nt = $V->sub_regions_sequence;
                $germ_aa = $V->sub_regions_translation;
                $self->_check_aa_nt_shared_regions($germ_nt,$germ_aa);
            }
        }
        # Check query V region;
        if ( $self->rendered_alignment && $self->rendered_alignment->query ) {
            my $V = $self->rendered_alignment->query;
            if ( $V->sub_regions_sequence && $V->sub_regions_translation ) {
                $query_nt = $V->sub_regions_sequence;
                $query_aa = $V->sub_regions_translation;
                $self->_check_aa_nt_shared_regions($query_nt,$query_aa);
            }
        }


        if ( $germ_nt && $germ_aa && $query_nt && $query_aa ) {
            if ( $self->is_almost_perfect ) {
                $mismatches = $self->_mutation_relative_to_germiline_nt($query_nt, $germ_nt, $query_aa, $germ_aa );
            }
        }
        return $mismatches;
    }


    method _check_aa_nt_shared_regions (Object $nt, Object $aa) {
        my @accessors_in_nt;
        my $meta = $nt->meta;
        # Get all accessor present in nt
        for my $attr ( $meta->get_all_attributes ) {
            my $accessor = $attr->name;
            my $predicate = "has_".$accessor;
            push @accessors_in_nt, $accessor if $nt->$predicate;
        }

        # Check if we have the same accessors in aa
        foreach my $accessor (@accessors_in_nt) {
            next if $accessor =~ /start|end/i;
            my $predicate = "has_".$accessor;
            unless ($aa->$predicate){
                print "query_id: ".$self->query_id."\n";
                print "Cannot find accessor in translation: ". $accessor."\n";
                p $nt;
                p $aa;
                die;
            }
        } 
    }


    method _mutation_relative_to_germiline_nt (Object $query_nt, Object $germ_nt, Object $query_aa, Object $germ_aa) {
        my @regions = (qw/ FWR1 CDR1 FWR2 CDR2 FWR3 /);
        
        my @germline_array = ('N') x ( $germ_nt->FWR1_start - 1 ); # keep diff compared to germline
        my @germline_array_aa = ('X') x (( $germ_nt->FWR1_start - 1 )/3);

        my @query_array = @germline_array; # keep query sequence
        my @query_array_aa = @germline_array_aa;

        my %hash_germline_regions;
        my %hash_query_regions;

        # aminoacid version
        my %hash_germline_regions_aa;
        my %hash_query_regions_aa;
        
        push @{$hash_germline_regions{FWR1}}, @germline_array;
        push @{$hash_query_regions{FWR1}}, @germline_array;
        
        foreach my $r (@regions) {
            my $predicate = 'has_' . $r;
            next unless $germ_nt->$predicate && $query_nt->$predicate;
            
            my ( $query, $germ ) = $self->_inspect_size( $query_nt->$r, $germ_nt->$r, $r );
            my ( $query_aa, $germ_aa ) = $self->_tranlaste( $query, $germ );

            my $germline_mutations = $self->_compare_string($query,$germ, use_dots => 1);
            my $query_sequence = $self->_compare_string($query,$germ, use_dots => 0 );
 
            my $germline_mutations_aa = $self->_compare_string($query_aa, $germ_aa, use_dots => 1);
            my $query_sequence_aa = $self->_compare_string($query_aa, $germ_aa, use_dots => 0 );
            
            push @germline_array, split '', $germline_mutations;
            push @query_array, split '', $query_sequence;
    
            push @germline_array_aa, split '', $germline_mutations_aa;
            push @query_array_aa, split '', $query_sequence_aa;
            
            push @{ $hash_germline_regions{$r} }, split '', $germline_mutations;
            push @{ $hash_query_regions{$r} },    split '', $query_sequence;

            push @{ $hash_germline_regions_aa{$r} }, split '', $germline_mutations_aa;
            push @{ $hash_query_regions_aa{$r} },    split '', $query_sequence_aa;
        }
        
        return {
            germ_regions => \%hash_germline_regions, 
            germ_regions_aa => \%hash_germline_regions_aa, 
            complete_germ => \@germline_array,
            complete_germ_aa => \@germline_array_aa,
        };
    }


    method _compare_string ($query, $germ, :$use_dots ) {
        my @g = split '', $germ;
        my @q = split '', $query;
        my $result;
        if ($use_dots){
            $result = join '', map { $g[$_] eq $q[$_] ? '.' : $q[$_] } 0 .. $#q;
        }
        else {
            $result = join '', map { $g[$_] eq $q[$_] ? $g[$_] : $q[$_] } 0 .. $#q;
        }
        return $result;
    }

    method _inspect_size ($query,$germ, $region) {
        my $l_germ  = length $germ;
        my $l_query = length $query;
        my $final_germ; 
        my $final_query;
        if ( $l_germ == $l_query ) {
            $final_germ = $germ;
            $final_query= $query;
        }
        else {
            if ( $l_query > $l_germ && $region =~ /FWR3/) {
                # Correct
                $final_germ = $germ;
                my @q = split '', $query;
                $final_query = join '', @q[0..($l_germ-1)]; 
            }
            else {
                p $germ;
                p $query;
                say $self->igblast_block;
                p $self;
                die "Cannot compare strings with different sizes"
 
            }
        }

        # remove insertion from germline
        if ( $final_germ =~ /\-/ ) {
            my @g = split '', $final_germ;
            my @q = split '', $final_query;
            my $aux_germ  = $g[0];
            my $aux_query = $q[0];
            for ( 1 .. $#q ) {
                $aux_query .= $q[$_] unless $g[$_] eq '-';
                $aux_germ  .= $g[$_] unless $g[$_] eq '-';
            }

            $final_germ  = $aux_germ;
            $final_query = $aux_query;
        }

        return ($final_query, $final_germ);
    }

    method _tranlaste ($query, $germ) {
        my $query_obj = Bio::PrimarySeq->new(
            -id => 'query',
            -seq => $query,
        );
        my $germ_obj = Bio::PrimarySeq->new(
            -id => 'germ',
            -seq => $germ,
        );
        
        return( $query_obj->translate->seq, $germ_obj->translate->seq );
    }

    # Test mismatches after object contruction;
    sub BUILD { my $self = shift; $self->mismatches; }

    method infer_CDR3_nt {
        my $cdr3_seq = "N/A";
    
        if ( $self->rendered_alignment ) {
            my $query = $self->rendered_alignment->query;

            if ( $query->sub_regions_sequence->has_FWR3 ) {

                my $r = quotemeta $query->sub_regions_sequence->FWR3;

                # If heavy chain
                if ( $self->chain_type =~ /heavy/i ) {

                    if ( $query->sequence =~ /($r)(\S+)(TGGGG[ATCG])/i ) {
                        $cdr3_seq = $2 . $3."|";
                    }
                    elsif ( $query->sequence =~ /($r)(\S+)/i ) {
                        $cdr3_seq = $2;
                    }
                }

                # If Iflight chain
                elsif ( $self->chain_type =~ /kappa|lambda/i ) {

                    if ( $self->database =~ /human/i ) {
                        if ( $query->sequence =~ /($r)(\S+)(TT[CT]GG[TCA])/i ) {
                            $cdr3_seq = $2 . $3."|";
                        }
                        elsif ( $query->sequence =~ /($r)(\S+)/i ) {
                            $cdr3_seq = $2;
                        }
                    }
                    elsif ( $self->database =~ /mouse/i ) {
                        if ( $query->sequence =~ /($r)(\S+)([ACGT][ACGT][ACGT]GG[TCA])[ACTG][ACGT][ACGT]GG[TGCA](AC[ACTG]AA[AG]){0,1}/i ) {
                            $cdr3_seq = $2 . $3. "|";
                        }
                        elsif ( $query->sequence =~ /($r)(\S+)/i ) {
                            $cdr3_seq = $2;
                        }
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

                if ( $self->chain_type =~ /heavy/i ) {
                    if ( $query->translation =~ /($r)(\S+)(WG)/i ) {
                        $cdr3_seq = $2 . $3. "|";
                    }
                    elsif ( $query->translation =~ /($r)(\S+)/i ) {
                        $cdr3_seq = $2;
                    }
                }

                # If light chain
                elsif ( $self->chain_type =~ /kappa|lambda/i ) {
                    if ( $self->database =~ /human/i ) {
                        if ( $query->translation =~ /($r)(\S+)(FG)/i ) {
                            $cdr3_seq = $2 . $3 . "|";
                        }
                        elsif ( $query->translation =~ /($r)(\S+)/i ) {
                            $cdr3_seq = $2;
                        }
                    }
                    elsif ( $self->database =~ /mouse/i ) {
                        if ( $query->translation =~ /($r)(\S+)(\S{1}G)\S{1}GTK/i ) {
                            $cdr3_seq = $2 . $3 . "|";
                        }
                        elsif ( $query->translation =~ /($r)(\S+)/i ) {
                            $cdr3_seq = $2;
                        }
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
                    $type = 'heavy';
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
                    elsif ( $self->chain_type =~ /heavy/i ) {
                        $answer = 1 if $r->top_D_match && $r->top_D_match ne 'N/A';
                    }
                }
            }
        }

        return $answer;
    }

    method is_v_complete {
        my $answer  = 0;
        my @regions = (qw/ FWR1 FWR2 FWR3 CDR1 CDR2 /);
        my @perfect_regions;
        if ( $self->rendered_alignment && $self->rendered_alignment->query ) {
            my $V = $self->rendered_alignment->query;
            if ( $V->sub_regions_sequence && $V->sub_regions_translation ) {
                my $query_nt = $V->sub_regions_sequence;
                foreach my $region (@regions) {
                    my $predicate = "has_" . $region;
                    if ( $query_nt->$predicate ) {
                        my $region_length = length $query_nt->$region;
                        my @count_n       = $query_nt->$region =~ m/N+/i;
                        my $ratio         = scalar @count_n / $region_length;
                        if ( $ratio > 0.2 ) {
                            return 0;
                        }
                        else {
                            push @perfect_regions, $region;
                        }
                    }
                    else {
                        return 0;
                    }
                }
            }

        }
        $answer = 1 if ( scalar @regions == scalar @perfect_regions );
        return $answer;
    }


    # almost perfect means is_complete and also it has all FRW and CDR regions
    # with 20% or less of nucleotides unknown
    method is_almost_perfect {
        my $answer      = 0;
        my @regions     = (qw/ FWR1 FWR2 FWR3 CDR1 CDR2 /);
        my @perfect_regions;
        if ( $self->is_complete ) {
            if ( $self->rendered_alignment && $self->rendered_alignment->query ) {
                my $V = $self->rendered_alignment->query;
                if ( $V->sub_regions_sequence && $V->sub_regions_translation ) {
                    my $query_nt = $V->sub_regions_sequence;
                    foreach my $region (@regions) {
                        my $predicate = "has_" . $region;
                        if ( $query_nt->$predicate ) {
                            my $region_length = length $query_nt->$region; 
                            my @count_n = $query_nt->$region =~ m/N+/i; 
                            my $ratio = scalar @count_n/$region_length;
                            if ( $ratio > 0.2 ){
                                return 0
                            }
                            else {
                                push @perfect_regions,$region;
                            }
                        }
                        else {
                            return 0;
                        }
                    }
                }
            }

        }
        $answer = 1 if (scalar @regions == scalar @perfect_regions);
        return $answer;
    }


    # perfect means is_complete and also it has all FRW and CDR regions + those regions are without N (good quality)
    method is_perfect {
        my $answer      = 0;
        my @regions     = (qw/ FWR1 FWR2 FWR3 CDR1 CDR2 /);
        my @perfect_regions;
        if ( $self->is_complete ) {
            if ( $self->rendered_alignment && $self->rendered_alignment->query ) {
                my $V = $self->rendered_alignment->query;
                if ( $V->sub_regions_sequence && $V->sub_regions_translation ) {
                    my $query_nt = $V->sub_regions_sequence;
                    foreach my $region (@regions) {
                        my $predicate = "has_" . $region;
                        if ( $query_nt->$predicate ) {
                            if ( $query_nt->$region =~ /N+/i ){
                                return 0
                            }
                            else {
                                push @perfect_regions,$region;
                            }
                        }
                        else {
                            return 0;
                        }
                    }
                }
            }

        }
        $answer = 1 if (scalar @regions == scalar @perfect_regions);
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
        if ($self->chain_type =~ /heavy/i){
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
