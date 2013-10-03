package TestsFor::Bio::Moose::IgBlast::Rearrangement;
    use Test::Class::Moose;
    use Bio::Moose::IgBlast::Rearrangement;

    sub test_construction {
        my ( $test, $report ) = @_;
        my $obj = Bio::Moose::IgBlast::Rearrangement->new(
            'top_V_match' => 'IGHV1-3*01',
            'top_D_match' => 'IGHD3-22*01,IGHD2-2*02,IGHD2-2*02',
            'top_J_match' => 'IGHV1-3*01',
            'chain_type'  => 'VH',
            'V_J_frame'   => 'In-frame',
            'strand'      => '-',
            'stop_codon'  => 'Yes',
            'productive'  => 'No',
        );
        isa_ok $obj, 'Bio::Moose::IgBlast::Rearrangement';

    }
1
