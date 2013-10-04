package TestsFor::Bio::Moose::IgBlast;
    use Test::Class::Moose;
    with 'Test::Class::Moose::Role::AutoUse';
    use Bio::Moose::IgBlast;
    
    BEGIN {
        has 'obj' => (
            is            => 'rw',
            isa           => 'Object',
            documentation => 'Hold object for test purposes',
        );
    }
    
    sub extra_constructor_args {}

    sub test_setup{
         my $test = shift;
         $test->obj($test->class_name->new(
            database              => 'filename',
            domain_classification => 'kabat',
            molecule              => 'N',
            query_id              => "id1",
            version               => '2.26+'
             )
         );
    }

    sub test_construction {
        my ( $test, $report ) = @_;
        isa_ok $test->obj, 'Bio::Moose::IgBlast';
    }

    sub test_accessors {
        my ( $test, $report ) = @_;
        my $o = $test->obj;
        is( $o->database,              'filename', "Read Database" );
        is( $o->domain_classification, 'kabat',    "Read Domain" );
        is( $o->molecule,              'N',        "Read molecule" );
        is( $o->query_id,              'id1',      "Read query_id" );
        is( $o->version,               '2.26+',    "Read version" );
    }
1
