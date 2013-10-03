use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::IgBlast {
    use MooseX::StrictConstructor;
    has 'molecule'              => ( is => 'ro', isa => 'Str',  required => 1 );
    has 'version'               => ( is => 'ro', isa => 'Str',  required => 1 );
    has 'query_id'              => ( is => 'ro', isa => 'Str',  required => 1 );
    has 'database'              => ( is => 'ro', isa => 'Str',  required => 1 );
    has 'domain_classification' => ( is => 'ro', isa => 'Str',  required => 0 );
    has 'converted_sequence'    => ( is => 'rw', isa => 'Bool', required => 1, default => 0 );
    has 'rearrangement_summary' => (
        is       => 'ro',
        isa      => 'Bio::Moose::IgBlast::Rearrangement',
        required => 0
    );
    has 'junction_details' => (
        is       => 'ro',
        isa      => 'Bio::Moose::IgBlast::Junction',
        required => 0
    );
    has 'alignments' => (
        is       => 'ro',
        isa      => 'Bio::Moose::IgBlast::Alignment',
        required => 0
    );
    has 'hit_table' => (
        is       => 'ro',
        isa      => 'Bio::Moose::IgBlast::HitTable',
        required => 0
    );

    # store any information you want in this part
    has misc     => ( is => 'rw', isa => 'Any' );
    has init_pos => ( is => 'ro', isa => 'Int' );
}

 # ABSTRACT: turns baubles into trinkets
