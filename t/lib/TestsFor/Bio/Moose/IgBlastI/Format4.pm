package TestsFor::Bio::Moose::IgBlastI::Format4;
    use Test::Class::Moose;
    use Bio::Moose::IgBlastI::Format4;

    sub test_construction {
        my ( $test, $report ) = @_;
        my $obj = Bio::Moose::IgBlastI::Format4->new(
            'file'  => '',
        );
        isa_ok $obj, 'Bio::Moose::IgBlastI::Format4';
    }

    sub test_construction {
        my ( $test, $report ) = @_;
        can_ok 'Bio::Moose::IgBlastI::Format4', 'new';
        my $obj = Bio::Moose::IgBlastI::Format4->new( 'file' => '', );
        isa_ok $obj, 'Bio::Moose::IgBlastI::Format4';
    }

    sub test_parse_version_1_0 {
        my ( $test, $report ) = @_;
        my $obj = Bio::Moose::IgBlastI::Format4->new(
            'file'  => 't/data/format4/input_files/1.0/IgBlast_heavy_lotta_1.0-4_all.out',
        );
        can_ok $obj,(qw/features next_feature/);
        isa_ok $obj->features, 'ARRAY';
        my $first = $obj->next_feature;
        isa_ok $first, 'Bio::Moose::IgBlast';
    }

    sub test_parse_version_1_1 {
        my ( $test, $report ) = @_;
        my $obj = Bio::Moose::IgBlastI::Format4->new(
            'file'  => 't/data/format4/input_files/1.1/IgBlast_heavy_lotta_1.1-4_all.out',
        );
        can_ok $obj,(qw/features next_feature/);
        isa_ok $obj->features, 'ARRAY';
        my $first = $obj->next_feature;
        isa_ok $first, 'Bio::Moose::IgBlast';
    }

    sub test_parse_version_1_2 {
        my ( $test, $report ) = @_;
        my $obj = Bio::Moose::IgBlastI::Format4->new(
            'file'  => 't/data/format4/input_files/1.2/IgBlast_heavy_lotta_1.2-4_all.out',
        );
        can_ok $obj,(qw/features next_feature/);
        isa_ok $obj->features, 'ARRAY';
        my $first = $obj->next_feature;
        isa_ok $first, 'Bio::Moose::IgBlast';
    }

1
