use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::IgBlast::Alignment::Region {
    use MooseX::StrictConstructor;
    has 'from'              => ( is => 'ro', isa => 'Int|Str', required => 1 );
    has 'to'                => ( is => 'ro', isa => 'Int|Str', required => 1 );
    has 'length'            => ( is => 'ro', isa => 'Int|Str', required => 1 );
    has 'matches'           => ( is => 'ro', isa => 'Int|Str', required => 1 );
    has 'mismatches'        => ( is => 'ro', isa => 'Int|Str', required => 1 );
    has 'gaps'              => ( is => 'ro', isa => 'Int|Str', required => 1 );
    has 'identinty_percent' => ( is => 'ro', isa => 'Int|Str', required => 1 );
}
