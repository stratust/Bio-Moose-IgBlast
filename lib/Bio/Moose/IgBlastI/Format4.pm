use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::IgBlastI::Format4 {
    use MooseX::StrictConstructor;
    use Bio::Moose::IgBlast;
    use Bio::Moose::IgBlast::Rearrangement;
    use Bio::Moose::IgBlast::Junction;
    use Bio::Moose::IgBlast::Alignment;
    use Bio::Moose::IgBlast::Alignment::Region;
    use Bio::Moose::IgBlast::HitTable;
    use Bio::Moose::IgBlast::HitTable::Entry;
    use Data::Printer;
    
    has 'file' => (
        is            => 'ro',
        isa           => 'Any',
        required      => 1,
        documentation => 'IgBlast file to be open',
    );
    
    has 'blast_version' => (
        is            => 'rw',
        isa           => 'Str',
    );

    has 'features' => (
        is            => 'rw',
        isa           => 'ArrayRef[Bio::Moose::IgBlast]',
        traits        => ['Array'],
        lazy          => 1,
        builder       => '_build_features',
        documentation => 'ArrayRef of features',
        handles       => {
            all_features   => 'elements',
            add_features   => 'push',
            next_feature   => 'shift',
            map_features   => 'map',
            count_features => 'count',
        },
    );


    method _build_features {
        my @objects;
        my $init_pos = 1;

        # By default Perl pulls in chunks of text up to a newline (\n) character; newline is
        # the default Input Record Separator. You can change the Input Record Separator by
        # using the special variable "$/". When dealing with IgBlast files I normally change the
        # Input Record Separator to "IGBLASTN" which allows your script to take in a full, multiline
        # IGBLAST record at once
        $/ = "IGBLAST";

        open( my $in, '<', $self->file ) || die "Cannot open/read file " . $self->file . "!";

        my $junk = <$in>;    # Discard the first "IGBLAST"

        while ( my $block = <$in> ) {
            # Remove "IGBLAST"
            chomp $block;

            # Each IgBLAST aligment (block) can be divide em sub_blocks:
            my ( $info_block, $rearrangement_block, $junction_block, $aln_block, $hit_table_block, $converted_block )
                = $self->_guess_sub_blocks($block);

            my %obj_params;

            # Each sub block has specific parser methods
            my ( $info_param, $rearrangement_param, $junction_param, $aln_param, $hit_table_param );
            if ($info_block) {
                $info_param = $self->_parse_info_block($info_block);
                %obj_params = (%obj_params, %{$info_param});
            }
            if ($rearrangement_block) {
                $rearrangement_param = $self->_parse_rearrangement_block($rearrangement_block);
                %obj_params = (%obj_params, %{$rearrangement_param});
            }
            if ($junction_block) {
                $junction_param = $self->_parse_junction_block($junction_block);
                %obj_params = (%obj_params, %{$junction_param});
            }
            if ($aln_block) {
                $aln_param = $self->_parse_alignment_block($aln_block);
                %obj_params = (%obj_params, %{$aln_param});
            }
            if ($hit_table_block) {
                $hit_table_param = $self->_parse_hit_table_block($hit_table_block);
                %obj_params = (%obj_params, %{$hit_table_param});
            }
    

            my $obj =
                Bio::Moose::IgBlast->new( %obj_params, init_pos => $init_pos );

            if ($converted_block){
                $obj->converted_sequence(1);
            }

            push( @objects, $obj );
            $init_pos++;
        }
        close($in);
        
        return \@objects;
    }


    method _guess_sub_blocks (Str $block) {
        # IGBLAST record at once
        local $/ = "IGBLAST";
        chomp $block;

        # Each IgBLAST aligment (block) can be divide em sub_blocks:
        my ( $info_block, $rearrangement_block, $junction_block, $aln_block, $hit_table_block , $converted_block);
        
        my @sub_blocks = split /^\s*$/m, $block;

        while (my $sb = shift @sub_blocks){
            if ( $sb =~ /Query/mi && $sb =~ /domain classification/mi ) {
                $info_block = $sb;
            }
            elsif ( $sb =~ /rearrangement summary/mi ) {
                $rearrangement_block = $sb;
            }
            elsif ( $sb =~ /junction details/mi ) {
                $junction_block = $sb;
            }
            elsif ( $sb =~ /alignment/mi && $sb =~ /top germline/ ) {
                $aln_block = $sb;
            }
            elsif ( $sb =~ /hit table/mi && $sb =~ /fields/mi ) {
                $hit_table_block = $sb;
            }
            elsif ( $sb =~ /converted to the plus strand/mi) {
                $converted_block = $sb;
            }
            else{
                die "Cannot guess sub-block type for: ".$sb;
            }
        } 

        return ( $info_block, $rearrangement_block, $junction_block, $aln_block, $hit_table_block, $converted_block );
    }


    method _parse_info_block ($info) {
        my %hash;
        if ($info =~ /^
            ([NP])          # get the molecule
            \s+
            (\S+[\+])       # get IgBlast version
            .*Query:\s* 
            (\S+)           # Get query name
            .*Database:\s* 
            (.*)\n            # Get database V D and J
            \W*\#\s+Domain\s+classification\s+requested:\s*
            (\S+)           # Get classification
            /xs
            )
        {
            %hash = (
                molecule              => $1,
                version               => $2,
                query_id              => $3,
                database              => $4,
                domain_classification => $5,
            );
            $self->blast_version($2);
        }
        else{
            die "Problem with info_block :\n$info";
        }
        return \%hash;
    }


    method _parse_rearrangement_block (Str $rearrangement) {
        my ( %hash, @fields );
        local $/ = "\n";
        open( my $in, '<', \$rearrangement );

        while ( my $row = <$in> ) {
            chomp $row;
            next if $row =~ /^$/;

            if ( $row =~ /^\#.*\((.*)\)/ ) {
                @fields = split /\s*,\s*/, $1;
                
                # normalize fields names
                foreach my $f (@fields) {
                    if ( $f =~ /Top V/i ) {
                        $f = 'top_V_match';
                    }
                    elsif ( $f =~ /Top J/i ) {
                        $f = 'top_J_match';
                    }
                    elsif ( $f =~ /Top D/i ) {
                        $f = 'top_D_match';
                    }
                    elsif ( $f =~ /chain type/i ) {
                        $f = 'chain_type';
                    }
                    elsif ( $f =~ /V\S{1}J\s+frame/i ) {
                        $f = 'V_J_frame';
                    }
                    elsif ( $f =~ /productive/i ) {
                        $f = 'productive';
                    }
                    elsif ( $f =~ /stop\s*codon/i ) {
                        $f = 'stop_codon';
                    }
                    elsif ( $f =~ /strand/i ) {
                        $f = 'strand';
                    }
                }
                next;
            }
            else {
                my @values = split /\t/, $row;
                if ( $#fields == $#values ) {
                    my %aux;
                    @aux{@fields} = @values;
                    $hash{rearrangement_summary} =
                        Bio::Moose::IgBlast::Rearrangement->new( %aux );
                }
                else {
                    p @fields;
                    p @values;
                    die "Problem with rearrangement_block :\n$rearrangement";

                }

            }

        }
        close($in);
        return \%hash;
    }


    method _parse_junction_block (Str $junction) {
        my ( %hash, @fields );
        local $/ = "\n";
        open( my $in, '<', \$junction );

        while ( my $row = <$in> ) {
            chomp $row;
            next if $row =~ /^$/;
            
            my $regex = qr/^[\#]\s*V.*[\(]D[\)].*J.*?\((.*)?\)\./;
            if ( $self->blast_version eq '2.2.26+' ) {
                $regex = qr/^[\#]\s*V.*[\(]D[\)].*J.*?\((.*J\s+start)\./;
            }
            if ( $row =~ /$regex/i ) {
                @fields = split /\s*,\s*/, $1;
                
                # normalize fields names
                foreach my $f (@fields) {
                    if ( $f =~ /V end/i ) {
                      $f = 'V_end';
                    }
                    elsif ( $f =~ /V\S{1}D\s+junction/i ) {
                        $f = 'V_D_junction';
                    }
                    elsif ( $f =~ /D\s+region/i ) {
                        $f = 'D_region';
                    }
                    elsif ( $f =~ /D\S{1}J\s+junction/i ) {
                        $f = 'D_J_junction';
                    }
                    elsif ( $f =~ /J\s+start/i ) {
                        $f = 'J_start';
                    }
                    elsif ( $f =~ /V\S{1}J\s+junction/i ) {
                        $f = 'V_J_junction';
                    }
                }
                next;
            }
            else {
                my @values = split /\s+/, $row;
                
                if ( $#fields == $#values ) {
                    my %aux;
                    @aux{@fields} = @values;
                    $hash{junction_details} =
                        Bio::Moose::IgBlast::Junction->new( %aux );
                }
                else {
                    die "Problem with junction_block :\n$junction";
                }

            }

        }
        close($in);
        return \%hash;

    }


    method _parse_alignment_block (Str $aln) {
        my %hash;
        local $/ = "\n";
        open( my $in, '<', \$aln );
        while ( my $row = <$in> ) {
            chomp $row;
            next if $row =~ /^$/;
            next if $row =~ /^\#/;
            $row =~ s/\s+\(.*\)\s+/ /g;
            my @f = split /\s+/, $row;

            die "Error with entry $row",p $aln if scalar @f != 8;

            my $region = Bio::Moose::IgBlast::Alignment::Region->new(
                from              => $f[1],
                to                => $f[2],
                length            => $f[3],
                matches           => $f[4],
                mismatches        => $f[5],
                gaps              => $f[6],
                identinty_percent => $f[7],
            );

            if ( $row =~ /^FW*R(\d+)/ ) {
                $hash{"FWR$1"} = $region;
            }
            elsif ( $row =~ /^CDR(\d+)/ ) {
                $hash{"CDR$1"} = $region;
            }
            elsif ( $row =~ /^Total/ ) {
                $hash{"Total"} = $region;
            }

        }
        close($in);
        my $aln_obj = Bio::Moose::IgBlast::Alignment->new(%hash);
        return {alignments => $aln_obj};
    }


    method _parse_hit_table_block (Str $hit_table) {
        my (%hash, @fields);
        local $/ = "\n";
        open( my $in, '<', \$hit_table );
        while ( my $row = <$in> ) {
            chomp $row;
            next if $row =~ /^$/;
            next if $row =~ /^\#\s*$/;
            next if $row =~ /^[\#]\s*BLAST/i;
            next if $row =~ /hit table/i;
            next if $row =~ /hits/i;

            my @values = (
                qw/qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop/
            );

            my @keys = (
                'query id',
                'query gi',
                'query acc.',
                'query acc.ver',
                'query length',
                'subject id',
                'subject ids',
                'subject gi',
                'subject gis',
                'subject acc.',
                'subject acc.ver',
                'subject accs.',
                'subject length',
                'q. start',
                'q. end',
                's. start',
                's. end',
                'query seq',
                'subject seq',
                'evalue',
                'bit score',
                'score',
                'alignment length',
                '% identity',
                'identical',
                'mismatches',
                'positives',
                'gap opens',
                'gaps',
                '% positives',
                'query/sbjct frames',
                'query frame',
                'sbjct frame',
                'BTOP'
            );
            
            my %field_index;
            @field_index{@keys} = @values;
            
            if ( $row =~ /^[\#]\s*Fields\s*:\s*(.*)/i ) {
                @fields = split /\s*,\s*/, $1;

                # normalize fields names
                foreach my $f (@fields) {
                    $f = $field_index{$f} if $field_index{$f};
                }
                next;
            }
            else {
                my @values = split /\s+/, $row;
                # Fields don't include the region
                if ( ( $#fields + 1 ) == $#values ) {
                    my %aux;
                    my $region = shift @values;
                    @aux{@fields} = @values;
                    my $entry = Bio::Moose::IgBlast::HitTable::Entry->new( %aux );
                    push @{ $hash{$region} }, $entry;
                }
                else {
                    p $row;
                    p @fields;
                    p @values;
                    die "Problem with hit_table_block :\n$hit_table";
                }

            }

        }

        my $hit_obj = Bio::Moose::IgBlast::HitTable->new(%hash);
        return { hit_table => $hit_obj };
    }

}
