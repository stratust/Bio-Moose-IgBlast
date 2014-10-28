use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::IgBlast::GermlineScore {
    use MooseX::StrictConstructor;
    has 'sequence_name' => ( is => 'ro', isa => 'Str', required => 1 );
    has 'score'         => ( is => 'ro', isa => 'Num', required => 1 );
    has 'evalue'        => ( is => 'ro', isa => 'Num', required => 1 );
}
