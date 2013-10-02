use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::IgBlast::Junction {
    use MooseX::StrictConstructor;
    has 'V_end'        => ( is => 'ro', isa => 'Str', required => 1 );
    has 'V_D_junction' => ( is => 'ro', isa => 'Str', required => 0 );
    has 'D_region'     => ( is => 'ro', isa => 'Str', required => 0 );
    has 'D_J_junction' => ( is => 'ro', isa => 'Str', required => 0 );
    has 'J_start'      => ( is => 'ro', isa => 'Str', required => 1 );
    # Kappa chain only
    has 'V_J_junction' => ( is => 'ro', isa => 'Str', required => 0 );
}
