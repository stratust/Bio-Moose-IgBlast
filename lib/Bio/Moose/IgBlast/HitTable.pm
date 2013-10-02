use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::IgBlast::HitTable {
    use MooseX::StrictConstructor;
    use MooseX::Aliases;
    has 'V' => (
        is       => 'ro',
        isa      => 'ArrayRef[Bio::Moose::IgBlast::HitTable::Entry]',
        required => 1,
        alias    => [qw/VK VH/],
    );
    has 'D' => (
        is       => 'ro',
        isa      => 'ArrayRef[Bio::Moose::IgBlast::HitTable::Entry]',
        required => 0,
        alias    => [qw/DK DH/],
    );
    has 'J' => (
        is       => 'ro',
        isa      => 'ArrayRef[Bio::Moose::IgBlast::HitTable::Entry]',
        required => 0,
        alias    => [qw/JK JH JL/],
    );
}
