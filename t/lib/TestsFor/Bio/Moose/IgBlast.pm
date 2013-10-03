package TestsFor::Bio::Moose::IgBlast;
    use Test::Class::Moose;
    use Bio::Moose::IgBlast;

    sub test_construction {
        my ( $test, $report ) = @_;
        my $obj = Bio::Moose::IgBlast->new(
            database              => 'filename',
            domain_classification => 'kabat',
            molecule              => 'N',
            query_id              => "id1",
            version               => '1'
        );
        isa_ok $obj, 'Bio::Moose::IgBlast';
    }
1
