use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::IgBlast::Alignment {
    use MooseX::StrictConstructor;
    has 'FWR1'  => ( is => 'ro', isa => 'Bio::Moose::IgBlast::Alignment::Region', required => 0 );
    has 'CDR1'  => ( is => 'ro', isa => 'Bio::Moose::IgBlast::Alignment::Region', required => 0 );
    has 'FWR2'  => ( is => 'ro', isa => 'Bio::Moose::IgBlast::Alignment::Region', required => 0 );
    has 'CDR2'  => ( is => 'ro', isa => 'Bio::Moose::IgBlast::Alignment::Region', required => 0 );
    has 'FWR3'  => ( is => 'ro', isa => 'Bio::Moose::IgBlast::Alignment::Region', required => 0 );
    has 'CDR3'  => ( is => 'ro', isa => 'Bio::Moose::IgBlast::Alignment::Region', required => 0 );
    has 'Total' => ( is => 'ro', isa => 'Bio::Moose::IgBlast::Alignment::Region', required => 1 );
}

