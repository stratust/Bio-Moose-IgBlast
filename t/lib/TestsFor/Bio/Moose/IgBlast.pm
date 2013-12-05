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
        is( $o->infer_CDR3_nt,         'N/A',      "Infer CDR3 nt" );
        is( $o->infer_CDR3_nt_length,  'N/A',      "Infer CDR3 nt lenght" );
        is( $o->infer_CDR3_aa,         'N/A',      "Infer CDR3 aa" );
        is( $o->infer_CDR3_aa_length,  'N/A',      "Infer CDR3 aa lenght" );
    }
1
