use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::IgBlast::Rearrangement {
    use MooseX::StrictConstructor;
    has 'top_V_match' => ( is => 'ro', isa => 'Str',  required => 1 );
    has 'top_D_match' => ( is => 'ro', isa => 'Str',  required => 0 );
    has 'top_J_match' => ( is => 'ro', isa => 'Str',  required => 1 );
    has 'chain_type'  => ( is => 'ro', isa => 'Str',  required => 1 );
    has 'V_J_frame'   => ( is => 'ro', isa => 'Str',  required => 1 );
    has 'strand'      => ( is => 'ro', isa => 'Str',  required => 1 );
    has 'stop_codon'  => ( is => 'ro', isa => 'Str',  required => 0 );
    has 'productive'  => ( is => 'ro', isa => 'Str',  required => 0 );
}
