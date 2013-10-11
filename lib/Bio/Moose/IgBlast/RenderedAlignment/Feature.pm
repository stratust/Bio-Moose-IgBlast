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
}
