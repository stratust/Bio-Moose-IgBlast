use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::IgBlast::RenderedAlignment::Feature::Region {
    use MooseX::StrictConstructor;
    has 'FWR1'       => ( is => 'ro', isa => 'Str', required => 0, lazy_build => 1 );
    has 'FWR1_start' => ( is => 'ro', isa => 'Int', required => 0, lazy_build => 1 );
    has 'FWR1_end'   => ( is => 'ro', isa => 'Int', required => 0, lazy_build => 1 );
    
    has 'FWR2'       => ( is => 'ro', isa => 'Str', required => 0, lazy_build => 1 );
    has 'FWR2_start' => ( is => 'ro', isa => 'Int', required => 0, lazy_build => 1 );
    has 'FWR2_end'   => ( is => 'ro', isa => 'Int', required => 0, lazy_build => 1 );
 
    has 'FWR3'       => ( is => 'ro', isa => 'Str', required => 0, lazy_build => 1 );
    has 'FWR3_start' => ( is => 'ro', isa => 'Int', required => 0, lazy_build => 1 );
    has 'FWR3_end'   => ( is => 'ro', isa => 'Int', required => 0, lazy_build => 1 );

    has 'CDR1'       => ( is => 'ro', isa => 'Str', required => 0, lazy_build => 1 );
    has 'CDR1_start' => ( is => 'ro', isa => 'Int', required => 0, lazy_build => 1 );
    has 'CDR1_end'   => ( is => 'ro', isa => 'Int', required => 0, lazy_build => 1 );

    has 'CDR2'       => ( is => 'ro', isa => 'Str', required => 0, lazy_build => 1 );
    has 'CDR2_start' => ( is => 'ro', isa => 'Int', required => 0, lazy_build => 1 );
    has 'CDR2_end'   => ( is => 'ro', isa => 'Int', required => 0, lazy_build => 1 );
  
}
