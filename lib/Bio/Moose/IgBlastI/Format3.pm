use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::IgBlastI::Format3 {
    use MooseX::StrictConstructor;
    use Bio::Moose::IgBlast;
    use Bio::Moose::IgBlast::Rearrangement;
    use Bio::Moose::IgBlast::GermlineScore;
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

            # Each IgBLAST aligment (block) can be divide in sub_blocks:
            #
            my ( $info_block, $germline_block, $rearrangement_block, $junction_block, $aln_block, $rendered_aln_block , $converted_block)
                = $self->_guess_sub_blocks($block);

            my %obj_params;


            # Each sub block has specific parser methods
            my ( $info_param, $germline_param, $rearrangement_param, $junction_param, $aln_param, $rendered_aln_param );

            %obj_params = (%obj_params, %{$header_param});

            if ($info_block) {
                $info_param = $self->_parse_info_block($info_block);
                %obj_params = (%obj_params, %{$info_param});
            }
            if ($germline_block) {
                $germline_param = $self->_parse_germline_block($germline_block);
                %obj_params = (%obj_params, %{$germline_param});
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
                Bio::Moose::IgBlast->new( %obj_params, init_pos => $init_pos, igblast_block => $block );

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
        elsif ( $header =~ /^Database\:\s*(.*)\n\s*.*sequences/s ) {
            %hash = (
                database => $1,
                version => '1.4.0', # hard coded because there isnt info anymore
                molecule => 'N', # hardcoded to N
            );
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
            elsif ( $sb =~ /Length/ && $sb =~ /alignment/mi && $sb =~ /producing/mi ) {
                $info_block .= $sb;
                $germline_block = shift @sub_blocks;
            }
            elsif ( $sb =~ /Length/ && $sb !~ /alignment/mi && $sb !~ /producing/mi ) {
                $info_block .= $sb;
            }
            elsif ( $sb =~ /Length=/mi ) {
                $info_block .= $sb;
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
            \s+.*Length=\s*
            (\S+)
            .*
            \s+Domain\s+classification\s+requested:\s*
            (\S+)           # Get classification
            /xs
            )
        {
            %hash = (
                query_id              => $1,
                query_length          => $2,
                domain_classification => $3,
            );
        }
        
        elsif ($info =~ /^
            .*ThisQuery=\s* 
            (\S+)           # Get query name
            \s+.*Length=\s*
            (\S+)
            /xs
            )
        {
            %hash = (
                query_id              => $1,
                query_length          => $2,
            );
        }

        else{
            die "Problem with info_block :\n$info";
        }
        $query_id = $hash{query_id};
        return \%hash;
    }


    method _parse_germline_block (Str $germ) {
        my %hash;
        local $/ = "\n";
        open( my $in, '<', \$germ );
        my @aux;
        while ( my $row = <$in> ) {
            chomp $row;
            next if $row =~ /^$/;
            next if $row =~ /^\#/;
            my ($seq, $score, $evalue);
            ($seq, $score, $evalue) = ($1, $2, $3) if ( $row =~ /^(.*)\s+(\d+\S+)\s+(\d+\S+)/ ) ;

            my $entry = Bio::Moose::IgBlast::GermlineScore->new(
                sequence_name => $seq,
                score         => $score,
                evalue        => $evalue,
            );
            push @aux, $entry;
        }

        close($in);
        return { germline_score => \@aux };
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
            if ( $self->blast_version && $self->blast_version eq '2.2.26+' ) {
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
            
            if (scalar @f != 8){
                #"Error with entry $row",p $aln;
            }

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
        my $regions_found = 0; # flag if there is regions alredy and keep colecting sequences in order to infer CDR3 later
        
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
            my $fixed_space;
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
                        $regions_found = 1;
                    }
                    else {
                        die "Problem with region: $row from $rendered_aln";
                    }
                }
                # skip block if no region annotation if found
                elsif ( $i == 0 &&  $row !~ /\%/ && $aux[1] !~ /[\-\<\>]/  && $regions_found == 0){
                    next BLOCKS;
                }
                # Query row
                elsif ( $aux[2] && ($aux[0] eq '' && $aux[2] =~ /\d+/  ))
                {
                    my ($this_seq,$start,$end);
                    my $l_space = length($space);
                    
                    # If V regions don't start from the beginig of the alignment,
                    # then use translation to get correct space
                    if ($fixed_space) {
                        $this_seq = $1 if $row =~ /^.{$fixed_space}(\S+)/;
                        if ($this_seq) {
                            my $nt = $this_seq;

                            # remove insertions
                            $nt =~ s/[^acgtn]//gi;
                            my $orig_length = $aux[4] - $aux[2] + 1;
                            $start = $aux[2] + ( $orig_length - length($nt) );
                        }

                        $end = $aux[4];
                    }
                    # otherwise use the space given by the region
                    elsif ( $l_space && $l_space > 1 ) {
                        $this_seq = $1 if $row =~ /^.{$l_space}(\S+)/;
                        if ($this_seq) {
                            my $nt = $this_seq;
                            
                            # remove insertions
                            $nt =~ s/[^acgtn]//gi;
                            my $orig_length = $aux[4] - $aux[2] + 1;
                            $start = $aux[2] + ( $orig_length - length($nt) );
                        }

                    }
                    else {
                        $this_seq = $aux[3];
                        $start = $aux[2];
                        $end   = $aux[4];
                    }
                    unless ($this_seq){
                        next BLOCKS;
                        #print "$aln_block\n";
                        #print "$row\n";
                        #p @aux;
                        #print $this_seq;
                        #die "Cannot find sequence for query: " .$query_id;
                    }
                    $seq{query}{id} = $aux[1];
                #    push( @{ $seq{query}{starts} }, $start );
                    #push( @{ $seq{query}{ends} },   $end );
                    #push @{ $seq{query}{seq} }, $seq;
                    $seq{query}{starts}{$b_count} = $start;
                    $seq{query}{ends}{$b_count} = $end;
                    $seq{query}{seq}{$b_count} = $this_seq;
                }

                # germline row
                elsif ( $aux[0] =~ /[VDJ]/ && $aux[1] =~ /\%/ && $space) {                    
                    my ($this_seq,$start,$end);

                    my $l_space = length($space);

                    if ($fixed_space) {
                        $this_seq = $1 if $row =~ /^.{$fixed_space}(\S+)/;
                        if ($this_seq){
                            my $nt = $this_seq;
                            # remove insertions
                            $nt =~ s/[^acgtn]//gi;
                            my $orig_length = $aux[6] - $aux[4] + 1;
                            $start = $aux[4] + ($orig_length - length($nt));
                        }
                        $end = $aux[6];
                    }
                    # otherwise use the space given by the region
                    elsif ( $l_space > 1 ) {
                        $this_seq = $1 if $row =~ /^.{$l_space}(\S+)/;
                        if ($this_seq){
                            my $nt = $this_seq;
                            # remove insertions
                            $nt =~ s/[^acgtn]//gi;
                            my $orig_length = $aux[6] - $aux[4] + 1;
                            $start = $aux[4] + ($orig_length - length($nt));
                        }
                        $end = $aux[6];
                    }

                    else {
                        $this_seq   = $aux[5];
                        $start = $aux[4];
                        $end   = $aux[6];
                    }
                    die "Cannot find sequence" unless $this_seq;
                    
                     my $germ_id = $aux[3];
                     my $germ_type = $aux[0];
                     # sometimes igblast choose the same sequence (e.g IGHV5)
                     # for first and second best alignment so the to should be a
                     # combination of 
                     my $uniq_germ_id = "$germ_id|$aux[2]|$aux[1]";
 
                    # keep ids of best V D and J germlines
                    $seq{top_germline}{$germ_type} = $uniq_germ_id unless $seq{top_germline}{$germ_type};

                    $seq{germline}{$uniq_germ_id}{type}             = $germ_type;
                    $seq{germline}{$uniq_germ_id}{percent}          = $aux[1];
                    $seq{germline}{$uniq_germ_id}{identity}         = $aux[2];
                    $seq{germline}{$uniq_germ_id}{starts}{$b_count} = $start;
                    $seq{germline}{$uniq_germ_id}{ends}{$b_count}   = $end;
                    $seq{germline}{$uniq_germ_id}{seq}{$b_count}    = $this_seq;
                }

                # Translation
                else {
                    my $l_space = length $space;
                    my $aa;
                    if ($l_space) {
                        $aa = $1 if $row =~ /^\s{$l_space}(.*)/;
                        # Try harder to get the right translation
                        # when the regions don't start in the beginning
                        # of alignment
                        unless ($aa){
                            $aa = $1 if $row =~ /^.{$l_space}(.*)/;
                            $fixed_space = $l_space if $aa;
                        }

                        if ($second_translation) {
                            $seq{germline}{translation} .= $aa if $aa;
                        }
                        else {
                            $seq{query}{translation} .= $aa if $aa;
                        }
                        $second_translation = 1;
                    }
                    else {
                        # sometimes igblast show an complete alignment without
                        # a region defined, so is better skip this alignment
                        # block it anyway
                        next ROWS;
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
            print "Sequence_id: ".$query_id."\n"  unless $seq{query}{translation};
            my $q_subregion_seq = $self->_get_region_sequences( 
                $seq{regions}, 
                $seq{query}{seq},
                $seq{query}{starts},
                translation => $seq{query}{translation} );
            
            $aux_hash{query} =
                Bio::Moose::IgBlast::RenderedAlignment::Feature->new( %{$q_subregion_seq},
                id => $seq{query}{id}, );

            # For germline
            # V region
            unless ( $seq{germline}{translation} ) {
                print "Sequence_id: " . $query_id . "\n";
                p %seq;
            }

            my $germ_subregion_seq = $self->_get_region_sequences(
                $seq{regions},
                $seq{germline}{$best_V_id}{seq},
                $seq{germline}{$best_V_id}{starts},
                translation => $seq{germline}{translation}
            );

            $best_V_id =~ s/\|.*//g; 
            $aux_hash{best_V} =
                Bio::Moose::IgBlast::RenderedAlignment::Feature->new( %{$germ_subregion_seq},
                id => $best_V_id, );
             
           
            # D region
           my $best_D_id = $seq{top_germline}{D};

           if ($best_D_id) {
               $best_D_id =~ s/\|.*//g; 
 
                $aux_hash{best_D} =
                    Bio::Moose::IgBlast::RenderedAlignment::Feature->new( 
                    id => $best_D_id, );
            }

            # J region
            my $best_J_id = $seq{top_germline}{J};
   
            if ($best_J_id) {
                $best_J_id =~ s/\|.*//g; 
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


    method _get_region_sequences (Str $v_regions_str, ArrayRef | HashRef $seq_ref, HashRef $starts, Str : $translation) {
        my $seq_str_orig;
        my $seq_str;
        
        if ( ref $seq_ref eq 'ARRAY' ) {
            $seq_str = join '', @{$seq_ref};
        }
        else {
            foreach my $k ( sort { $a <=> $b } keys %{$seq_ref} ) {
                $seq_str .= $seq_ref->{$k};
            }
        }

        $seq_str_orig = $seq_str;

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
        my @parts_nt;
        @parts_nt = $seq_str =~ /$regex/ if $regex;
        $seq_str =~ s/\s+//g;
        $_ =~ s/\s+//g foreach @parts_nt;

        # Sometimes translation has less space at the end than DNA, so lets add
        # spaces at the end before split
        while (length($translation) < $regex_length){
            $translation .= ' ';
        }
        my @parts_aa;
        @parts_aa = $translation =~ /$regex/ if $regex;
        $translation =~ s/\s+//g;
        $_ =~ s/\s+//g foreach @parts_aa;

        my %seq_nt;
        my %seq_aa;
        @seq_nt{@r_name} = @parts_nt;
        @seq_aa{@r_name} = @parts_aa;
       
        # Correct start
        my @aux_starts = sort {$a <=> $b} values %{$starts};
        my $real_start = $aux_starts[0];

        my $i=0;
        my $removed = 0;
        my %seq_nt_range;
        foreach my $part (@parts_nt) {
            next unless $part;
            if ( $seq_str_orig =~ m/$part?/ ) {
                
                my ($region_start, $region_end);
                my $match_start = $-[0];
                my $match_end = $+[0];
                my $length = length($seq_str_orig) - $match_end;
                my $rest_of_sequence = substr $seq_str_orig, $match_end, $length;

                # find is are insertions in part
                my @insertion = $part =~ /\-/g;
                # for initial number insertion were already corrected;
                @insertion = () if $i == 0;
                
                # Igblast output is NOT ZERO-based!
                $region_start = $removed + $real_start + $match_start - scalar(@insertion);

                # Always remove 1 from END;
                $region_end = ( $removed + $real_start + $match_end ) - 1;

                #print $r_name[$i]."\n";
                #print $part."\n";
                #print "start: $region_start\n";
                #print "end: $region_end\n";
                
                $seq_nt_range{$r_name[$i]."_start"} = $region_start;
                $seq_nt_range{$r_name[$i]."_end"} = $region_end;

                
                # Keep the amount of sequence removed from original sequence;
                $removed += (length($seq_str_orig) - length($rest_of_sequence)) - scalar(@insertion);
                
                # Always use the rest of the sequence for next match (avoid
                # repetive regions always match the same place)
                $seq_str_orig = $rest_of_sequence;
            }
            $i++;
        }
       
        my $new_translation;
        foreach my $part_aa (@parts_aa) {
            $new_translation .= $part_aa if $part_aa;
        }
        my $new_seq;
        foreach my $part_nt (@parts_nt) {
            $new_seq .= $part_nt if $part_nt;
        }


        return {
            sub_regions_translation =>
                Bio::Moose::IgBlast::RenderedAlignment::Feature::Region->new(%seq_aa),
                sub_regions_sequence =>
                Bio::Moose::IgBlast::RenderedAlignment::Feature::Region->new(
                %seq_nt, %seq_nt_range
                ),
                translation         => $translation,
                translation_trimmed => $new_translation,
                sequence            => $seq_str,
                sequence_trimmed    => $new_seq
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
