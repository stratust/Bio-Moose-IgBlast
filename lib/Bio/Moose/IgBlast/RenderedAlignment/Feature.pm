use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::IgBlast::RenderedAlignment::Feature {
    use MooseX::StrictConstructor;
    has 'id'                      => ( is => 'ro', isa => 'Str',     required => 1 );
#    has 'identity'                => ( is => 'ro', isa => 'Str',     required => 1 );
#    has 'identity_percent'        => ( is => 'ro', isa => 'Str',     required => 1 );
    has 'translation'             => ( is => 'ro', isa => 'Str',     required => 0 );
    has 'sequence'                => ( is => 'ro', isa => 'Str',     required => 0 );
    has 'sub_regions_sequence'    => ( is => 'ro', isa => 'Bio::Moose::IgBlast::RenderedAlignment::Feature::Region', required => 0 );
    has 'sub_regions_translation' => ( is => 'ro', isa => 'Bio::Moose::IgBlast::RenderedAlignment::Feature::Region', required => 0 );

    method infer_CDR3_nt {
        my $cdr3_seq = "N/A";
        if ($self->sub_regions_sequence->has_FWR3){
            #print $self->sequence."\n";
            #die $self->sub_regions_sequence->FWR3;
            my $r = quotemeta $self->sub_regions_sequence->FWR3;
            if ($self->sequence =~ /($r)(\S+TGGGGC)/i){
                $cdr3_seq = $2."|";
            }
            elsif ($self->sequence =~ /($r)(\S+)/i){
                $cdr3_seq = $2;
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
        if ($self->sub_regions_translation->has_FWR3){
            #print $self->sequence."\n";
            #die $self->sub_regions_sequence->FWR3;
            my $r = quotemeta $self->sub_regions_translation->FWR3;
            if ($self->translation =~ /($r)(\S+WG)/i){
                $cdr3_seq = $2."|";
            }
            elsif ($self->translation =~ /($r)(\S+)/i){
                $cdr3_seq = $2;
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
