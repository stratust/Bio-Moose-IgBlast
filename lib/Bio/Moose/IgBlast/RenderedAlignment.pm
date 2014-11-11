use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::IgBlast::RenderedAlignment {
    use MooseX::StrictConstructor;
    use Data::Printer;

    has 'query' => (
        is            => 'rw',
        isa           => 'Bio::Moose::IgBlast::RenderedAlignment::Feature',
        required      => 1,
        lazy_build    => 1,
        documentation => 'Best V',
    );

    has 'best_V' => (
        is            => 'rw',
        isa           => 'Bio::Moose::IgBlast::RenderedAlignment::Feature',
        required      => 1,
        lazy_build    => 1,
        documentation => 'Best V',
    );

    has 'best_D' => (
        is            => 'rw',
        isa           => 'Bio::Moose::IgBlast::RenderedAlignment::Feature',
        required      => 0,
        lazy_build    => 1,
        documentation => 'Best D',
    );

    has 'best_J' => (
        is            => 'rw',
        isa           => 'Bio::Moose::IgBlast::RenderedAlignment::Feature',
        required      => 0,
        lazy_build    => 1,
        documentation => 'Best J',
    );

   
}
