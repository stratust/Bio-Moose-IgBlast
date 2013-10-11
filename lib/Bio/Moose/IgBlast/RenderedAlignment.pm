use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::IgBlast::RenderedAlignment {
    use MooseX::StrictConstructor;
    
    has 'query' => (
        is            => 'rw',
        isa           => 'Bio::Moose::IgBlast::RenderedAlignment::Feature',
        required      => 1,
        documentation => 'Best V',
    );  

    has 'best_V' => (
        is            => 'rw',
        isa           => 'Bio::Moose::IgBlast::RenderedAlignment::Feature',
        required      => 1,
        documentation => 'Best V',
    );

    has 'best_D' => (
        is            => 'rw',
        isa           => 'Bio::Moose::IgBlast::RenderedAlignment::Feature',
        required      => 0,
        documentation => 'Best D',
    );

    has 'best_J' => (
        is            => 'rw',
        isa           => 'Bio::Moose::IgBlast::RenderedAlignment::Feature',
        required      => 0,
        documentation => 'Best J',
    );
    
}
