use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::IgBlast::RenderedAlignment::Feature::Region {
    use MooseX::StrictConstructor;
    has 'FWR1' => ( is => 'ro', isa => 'Str', required => 0 );
    has 'FWR2' => ( is => 'ro', isa => 'Str', required => 0 );
    has 'FWR3' => ( is => 'ro', isa => 'Str', required => 0 );
    has 'CDR1' => ( is => 'ro', isa => 'Str', required => 0 );
    has 'CDR2' => ( is => 'ro', isa => 'Str', required => 0 );
}
