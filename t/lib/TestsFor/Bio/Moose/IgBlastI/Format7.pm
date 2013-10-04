package TestsFor::Bio::Moose::IgBlastI::Format7;
    use Test::Class::Moose;
    with 'Test::Class::Moose::Role::AutoUse';
    use Bio::Moose::IgBlastI::Format7;

    BEGIN {
        has 'obj' => (
            is            => 'rw',
            isa           => 'Object',
            documentation => 'Hold object for test purposes',
        );
    }

    sub test_setup {
        my $test = shift;
        $test->obj(
            $test->class_name->new(
                'file' =>
't/data/format7/input_files/1.0/IgBlast_heavy_lotta_1.0-7_all.out',
            )
        );
    }

    sub test_construction {
        my ( $test, $report ) = @_;

        isa_ok $test->obj, 'Bio::Moose::IgBlastI::Format7';
    }

    sub test_parse_version_1_0 {
        my ( $test, $report ) = @_;

        can_ok $test->obj,(qw/features next_feature/);
        
        isa_ok $test->obj->features, 'ARRAY';

        my $o = $test->obj->next_feature;
        isa_ok $o, 'Bio::Moose::IgBlast';

        # Test attributes
        is(
            $o->database,
            'database/human_gl_V database/human_gl_D database/human_gl_J',
            "Read Database"
        );
        is( $o->domain_classification, 'kabat', "Read Domain" );
        is( $o->molecule,              'N',     "Read molecule" );
        is( $o->query_id,              '3A_P94_LVH_G3-IgGinternal.ab1',   "Read query_id" );
        is( $o->version,               '2.2.27+', "Read version" );
        
        
        # Rearrangement Summary
        my $r = $o->rearrangement_summary;
        isa_ok $r, 'Bio::Moose::IgBlast::Rearrangement';

        # Junction Details
        my $j = $o->junction_details;
        isa_ok $j, 'Bio::Moose::IgBlast::Junction';

        # Alignments
        my $a = $o->alignments;
        isa_ok $a, 'Bio::Moose::IgBlast::Alignment';

        # Hit Table
        my $h = $o->hit_table;
        isa_ok $h, 'Bio::Moose::IgBlast::HitTable';

    }

1
