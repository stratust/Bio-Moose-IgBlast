use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::IgBlastI::Format3 {
    use MooseX::StrictConstructor;
    use Bio::Moose::IgBlast;
    use Bio::Moose::IgBlast::Rearrangement;
    use Bio::Moose::IgBlast::Junction;
    use Bio::Moose::IgBlast::Alignment;
    use Bio::Moose::IgBlast::Alignment::Region;
    use Bio::Moose::IgBlast::HitTable;
    use Bio::Moose::IgBlast::RenderedAlignment;
    use Bio::Moose::IgBlast::RenderedAlignment::Feature;
    use Bio::Moose::IgBlast::RenderedAlignment::Feature::Region;
    use Data::Printer;
    use Carp::Always;
    
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

    our $query_id;

    method _build_features {
        my @objects;
        my $init_pos = 1;

        # By default Perl pulls in chunks of text up to a newline (\n) character; newline is
        # the default Input Record Separator. You can change the Input Record Separator by
        # using the special variable "$/". When dealing with IgBlast files I normally change the
        # Input Record Separator to "IGBLASTN" which allows your script to take in a full, multiline
        # IGBLAST record at once
        local $/ = "Query=";

        open( my $in, '<', $self->file ) || die "Cannot open/read file " . $self->file . "!";

        my $header = <$in>;    # Discard the first "Query="
        chomp $header;
        my $header_param = $self->_parse_header($header);

        while ( my $block = <$in> ) {
            # Remove "Query"
            chomp $block;

            # Each IgBLAST aligment (block) can be divide em sub_blocks:
            #
            my ( $info_block, $germline_block, $rearrangement_block, $junction_block, $aln_block, $rendered_aln_block , $converted_block)
                = $self->_guess_sub_blocks($block);

            my %obj_params;

            # Each sub block has specific parser methods
            my ( $info_param, $rearrangement_param, $junction_param, $aln_param, $rendered_aln_param );

            %obj_params = (%obj_params, %{$header_param});

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

            if ($rendered_aln_block && $aln_block ) {
                $rendered_aln_param = $self->_parse_rendered_aln_block($rendered_aln_block);
                %obj_params = (%obj_params, %{$rendered_aln_param});
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


    method _parse_header ($header) {
        my %hash;
        local $/="\n";
        if ( $header =~ /^BLAST([NP])\s+(\S+).*Database\:\s*(.*)\n\s*.*sequences/s ) {
            %hash = (
                molecule =>$1,
                version => $2,
                database => $3,
            );
            $self->blast_version($2);
        }
        else{
            die "REGEX NOT WORKING";
        }
        return \%hash;
    }


    method _guess_sub_blocks (Str $block) {
        # IGBLAST record at once
        local $/ = "Query=";
        chomp $block;
        $block="ThisQuery=$block";
    
        # Each IgBLAST aligment (block) can be divide em sub_blocks:
        my ( $info_block, $germline_block, $rearrangement_block, $junction_block, $aln_block, $rendered_aln_block , $converted_block);
        
        my @sub_blocks = split /^\s*$/m, $block;
        
        while (my $sb = shift @sub_blocks){
            if ( $sb =~ /ThisQuery/mi ) {
                $info_block = $sb;
            }
            elsif ( $sb =~ /domain classification/mi ) {
                $info_block .= $sb;
            }
            elsif ( $sb =~ /rearrangement summary/mi ) {
                $rearrangement_block = $sb;
            }
            elsif ( $sb =~ /junction details/mi ) {
                $junction_block = $sb;
            }
            elsif ( $sb =~ /Length/ && $sb =~ /alignment/mi && $sb =~ /producing/mi ) {
                $germline_block = join "\n", ($sb, shift @sub_blocks);
            }
            elsif ( $sb =~ /Length/ && $sb !~ /alignment/mi && $sb !~ /producing/mi ) {
                $germline_block = $sb;
            }
            elsif ( $sb =~ /alignment/mi && $sb =~ /top germline/mi ) {
                $aln_block = $sb;
            }
            elsif ( $sb =~ /^Alignments\n/mi || $sb =~ /no hits found/mi ) {
                $rendered_aln_block = $sb;
                while ( my $rest = shift @sub_blocks ) {
                    $rendered_aln_block .= "\n$rest";
                }
            }
            elsif ( $sb =~ /converted to the plus strand/mi) {
                $converted_block = $sb;
            }
            else{
                die "Cannot guess sub-block type for: ".$sb;
            }
        } 

        return ( $info_block, $germline_block, $rearrangement_block, $junction_block, $aln_block, $rendered_aln_block , $converted_block);
    }


    method _parse_info_block ($info) {
        my %hash;
        if ($info =~ /^
            .*ThisQuery=\s* 
            (\S+)           # Get query name
            .*
            \s+Domain\s+classification\s+requested:\s*
            (\S+)           # Get classification
            /xs
            )
        {
            %hash = (
                query_id              => $1,
                domain_classification => $2,
            );
        }
        elsif ($info =~ /^
            .*ThisQuery=\s* 
            (\S+)           # Get query name
            /xs
            )
        {
            %hash = (
                query_id              => $1,
            );
        }

        else{
            die "Problem with info_block :\n$info";
        }
        $query_id = $hash{query_id};
        return \%hash;
    }


    method _parse_rearrangement_block (Str $rearrangement) {
        my ( %hash, @fields );
        local $/ = "\n";
        open( my $in, '<', \$rearrangement );

        while ( my $row = <$in> ) {
            chomp $row;
            next if $row =~ /^$/;

            if ( $row =~ /^V.*D.*\((.*)\)/ ) {
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
            
            my $regex = qr/V.*[\(]D[\)].*J.*?\((.*)?\)\./;
            if ( $self->blast_version eq '2.2.26+' ) {
                $regex = qr/V.*[\(]D[\)].*J.*?\((.*J\s+start)\./;
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
            next if $row =~ /Alignment summary/i;
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


    method _parse_rendered_aln_block ($rendered_aln) {
        my %hash;
        my @aln_blocks = split /\n\n/ms, $rendered_aln;
        local $/ = "\n";
        my %seq;
        my $b_count=0;
        my $space;

        BLOCKS: foreach my $aln_block (@aln_blocks) {

            next if $aln_block =~ /Alignments/;
            next if $aln_block =~ /Lambda/sm;
            next if $aln_block =~ /Matrix\:/sm;
            next if $aln_block =~ /Database\:/sm;
            next if $aln_block =~ /Effective/;

            open( my $in, '<', \$aln_block )
              || die "Cannot open/read file " . $aln_block . "!";

            # count block_row
            my $i         = 0;
            my $query_row = 'N/A';
            my $germ_pos = 0;
            my $second_translation = 0;
            ROWS: while ( my $row = <$in> ) {
                next if $row =~ /^$/;
                chomp $row;
                my @aux = split /\s+/, $row;
                if ( $row =~ /no hits found/ims ) {
                    last BLOCKS;
                }

                # First row in each block indicate regions;
                if ( $i == 0 && $row !~ /\%/ && $aux[1] =~ /[\-\<\>]/ ) {
                    if ( $row =~ /^(\s+)(\S+)/ ) {
                        my $region;
                        ( $space, $region ) = ( $1, $2 );
                        $seq{regions} .= $region;
                    }
                    else {
                        die "Problem with region: $row form $rendered_aln";
                    }
                }

                # Query row
                elsif ( $aux[2] && ($aux[0] eq '' && $aux[2] =~ /\d+/  ))
                {
                    $seq{query}{id} = $aux[1];
                    push( @{ $seq{query}{starts} }, $aux[2] );
                    push( @{ $seq{query}{ends} },   $aux[4] );
                    push @{ $seq{query}{seq} }, $aux[3];
                }

                # germline row
                elsif ( $aux[0] =~ /[VDJ]/ && $aux[1] =~ /\%/ ) {
                     my $germ_id = $aux[3];
                     my $germ_type = $aux[0];
 
                    # keep ids of best V D and J germlines
                    $seq{top_germline}{$germ_type} = $germ_id unless $seq{top_germline}{$germ_type};

                    $seq{germline}{$germ_id}{type}             = $germ_type;
                    $seq{germline}{$germ_id}{percent}          = $aux[1];
                    $seq{germline}{$germ_id}{identity}         = $aux[2];
                    $seq{germline}{$germ_id}{starts}{$b_count} = $aux[4];
                    $seq{germline}{$germ_id}{ends}{$b_count}   = $aux[6];
                    $seq{germline}{$germ_id}{seq}{$b_count}    = $aux[5];
                }

                # Translation
                else {
                    my $l_space = length $space;
                    my $aa;
                    if ($l_space) {
                        $aa = $1 if $row =~ /^\s{$l_space}(.*)/;
                        if ($second_translation) {
                            $seq{germline}{translation} .= $aa if $aa;
                        }
                        else {
                            $seq{query}{translation} .= $aa if $aa;
                        }
                        $second_translation = 1;
                    }
                    else {
                        print "no space found before regions row:" . p @aux;
                        die $rendered_aln;
                    }

                }

                $i++;
            }
            close($in);
            $b_count++;
        }
 
        my $best_V_id = $seq{top_germline}{V};
        my %aux_hash;
        if ($best_V_id) {

            # For query
            my $q_subregion_seq = $self->_get_region_sequences( $seq{regions}, $seq{query}{seq},
                translation => $seq{query}{translation} );

            $aux_hash{query} =
                Bio::Moose::IgBlast::RenderedAlignment::Feature->new( %{$q_subregion_seq},
                id => $seq{query}{id}, );

            # For germline
            # V region
            my $germ_subregion_seq = $self->_get_region_sequences(
                $seq{regions},
                $seq{germline}{$best_V_id}{seq},
                translation => $seq{germline}{translation}
            );

            $aux_hash{best_V} =
                Bio::Moose::IgBlast::RenderedAlignment::Feature->new( %{$germ_subregion_seq},
                id => $best_V_id, );
             
           
            # D region
            my $best_D_id = $seq{top_germline}{D};
            if ($best_D_id) {
                $aux_hash{best_D} =
                    Bio::Moose::IgBlast::RenderedAlignment::Feature->new( 
                    id => $best_D_id, );
            }

            # J region
            my $best_J_id = $seq{top_germline}{J};
            if ($best_J_id) {
                $aux_hash{best_J} =
                    Bio::Moose::IgBlast::RenderedAlignment::Feature->new( 
                    id => $best_J_id, );
            }

        }
        if (%aux_hash) {
            %hash = (
                rendered_alignment => Bio::Moose::IgBlast::RenderedAlignment->new( %aux_hash )
            );
        }
        return \%hash;
    }


    method _get_region_sequences (Str $v_regions_str, ArrayRef | HashRef $seq_ref, Str : $translation) {
        my $seq_str;
        if ( ref $seq_ref eq 'ARRAY' ) {
            $seq_str = join '', @{$seq_ref};
        }
        else {
            foreach my $k ( sort { $a <=> $b } keys %{$seq_ref} ) {
                $seq_str .= $seq_ref->{$k};
            }
        }

        $v_regions_str =~ s/^\<//;
        $v_regions_str =~ s/\>$//;

        # Correct beginning based on translation
        if ( $translation && $translation =~ /^(\s+)/ ) {
            my $space_size = length $1;
            if ( $space_size > 1 ) {
                my $rm = $space_size - 1;
                $v_regions_str = $1 if ( $v_regions_str =~ /^.{$rm}(.*)/ );
                $translation   = $1 if ( $translation =~ /^.{$rm}(.*)/ );
                $seq_str       = $1 if ( $seq_str =~ /^.{$rm}(.*)/ );
            }
        }
        
        my @regions = split '><', $v_regions_str;
        my @r_name;
        my $regex;
       
        my $regex_length = 0;
        foreach my $r (@regions) {

            # Getting region names and sizes;
            if ( $r =~ /([^-]+)/ ) {
                my $aux = $1;
                $aux =~ s/FR/FWR/g;
                if ( length($aux) == 4 ) {
                    push @r_name, $aux;
                    $regex .= '(.{' . ( length($r) + 2 ) . '})';
                }
                else {
                    $regex .= '.{' . ( length($r) + 2 ) . '}';
                }
                $regex_length += length($r) + 2;
            }
        }
        
        # For germline sometimes we have less sequence than regions, so let's
        # add spaces at the end before split
        while (length($seq_str) < $regex_length){
            $seq_str .= ' ';
        }
        my @parts_nt = $seq_str =~ /$regex/;
        $seq_str =~ s/\s+//g;
        $_ =~ s/\s+//g foreach @parts_nt;

        # Sometimes translation has less space at the end than DNA, so lets add
        # spaces at the end before split
        while (length($translation) < $regex_length){
            $translation .= ' ';
        }
        my @parts_aa = $translation =~ /$regex/;
        $translation =~ s/\s+//g;
        $_ =~ s/\s+//g foreach @parts_aa;

        my %seq_nt;
        my %seq_aa;
        @seq_nt{@r_name} = @parts_nt;
        @seq_aa{@r_name} = @parts_aa;

        return {
            sub_regions_translation =>
                Bio::Moose::IgBlast::RenderedAlignment::Feature::Region->new(%seq_aa),
            sub_regions_sequence =>
                Bio::Moose::IgBlast::RenderedAlignment::Feature::Region->new(%seq_nt),
            translation => $translation,
            sequence    => $seq_str
        };

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
