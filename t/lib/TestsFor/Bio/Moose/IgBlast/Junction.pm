package TestsFor::Bio::Moose::IgBlast::Junction;
    use Test::Class::Moose;
    use Bio::Moose::IgBlast::Junction;

    sub test_construction {
        my ( $test, $report ) = @_;
        my $obj = Bio::Moose::IgBlast::Junction->new(
            'V_end'        => '',
            'V_D_junction' => '',
            'D_region'     => '',
            'D_J_junction' => '',
            'J_start'      => '',
            'V_J_junction' => '',
        );
        isa_ok $obj, 'Bio::Moose::IgBlast::Junction';

    }
1
