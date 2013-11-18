#!/usr/bin/env perl
use Moose;
use feature qw(say);
use MooseX::Declare;
use Method::Signatures::Modifiers;
BEGIN { our $Log_Level = 'info' }
 
class MyApp is dirty {
    use MooseX::App qw(Color);
    use Log::Any::App '$log',
        -screen => { pattern_style => 'script_long' },
        -file => { path => 'logs/', level => 'debug' };

    has 'log' => (
        is            => 'ro',
        isa           => 'Object',
        required      => 1,
        default       => sub { return $log },
        documentation => 'Keep Log::Any::App object',
    );
}

class MyApp::Igblast2Excel {
    extends 'MyApp'; # inherit log
    use MooseX::App::Command;    # important
    use MooseX::FileAttribute;
    use lib '../lib'; 
    use Bio::Moose::IgBlastIO;
    use Data::Printer  deparse => 1, sort_keys => 1, class => {expand => 'all'};

    command_short_description q[This command is awesome];
    command_long_description q[This command is so awesome, yadda yadda yadda];

    has_file 'input_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(i)],
        must_exist    => 1,
        required      => 1,
        documentation => q[Very important option!],
    );

    method run {
        my $cmd;
        $cmd = $1 if __PACKAGE__ =~ /\:\:(.*)$/;
        $self->log->warn("==> Starting $cmd <==");

        # Code Here
        my $in = Bio::Moose::IgBlastIO->new(
            file   => $self->input_file,
            format => 'format4'
        );

        my $no_aln = 0;
        my $i      = 0;
        my @ROWS   = (qw/best_V best_D best_J/);
        my @regions = (qw/ FWR1 CDR1 FWR2 CDR2 FWR3 /);
        
        #my @col_header = (qw/query_id/, @ROWS ,qw/CDR3_seq CDR3_length FWR1 gaps mismatches CDR1 gaps mismatches FWR2 gaps mismatches CDR2 gaps mismatches FRW3 gaps mismatches/);
        my @col_header = (qw/query_id/,@ROWS ,qw/CDR3_nt CDR3_nt_length CDR3_aa CDR3_aa_length FWR1 gaps mismatches FWR1_aa CDR1 gaps mismatches CDR1_aa FWR2 gaps mismatches FWR2_aa CDR2 gaps mismatches CDR2_aa FRW3 gaps mismatches FRW3_aa/);
        say join "\t",@col_header;
        while ( my $feature = $in->next_feature ) {
            my @col;    #keep info to print
            my @region_col;

            #say "\n\n$i QUERY: " . $feature->query_id;
            push @col, $feature->query_id;

            # Getting Bio::Moose::IgBlast::RenderedAlignment object
            my $aln = $feature->rendered_alignment;

            if ($aln) {
                foreach my $row (@ROWS) {
                    my $pred = "has_$row";

                    if ( $aln->$pred ) {

                        #say "BEST $row:" . $aln->$row->id;
                        push @col, $aln->$row->id;

                        #my $sub_r = $aln->$row->sub_regions_sequence;
                        #my $sub_t = $aln->$row->sub_regions_translation;
                        #next unless $sub_r;

=cut 
                        if ( $row =~ /query/ ) {
                            foreach my $r (@regions) {
                                my $pred = "has_$r";
                                if ( $sub_r->$pred ) {
                                    say "$r: "
                                      . $sub_r->$r . " ("
                                      . length( $sub_r->$r ) . ")"
                                      . " [ mismatches $r: "
                                      . $feature->alignments->$r->mismatches
                                      . " ]";

                                    push @region_col, ( $sub_r->$r, $feature->alignments->$r->gaps, $feature->alignments->$r->mismatches );
                                }
                                say "$r: "
                                  . $sub_t->$r . " ("
                                  . length( $sub_t->$r ) . ")"
                                  if $sub_t->$pred;
                            }

                        }
=cut

                    }
                    else {
                        push @col, ("N/A");
                    }

                }
                my $sub_r = $aln->query->sub_regions_sequence;
                my $sub_t = $aln->query->sub_regions_translation;

                # CDR3
                push @col,(
                    $aln->query->infer_CDR3_nt,
                    $aln->query->infer_CDR3_nt_length,
                    $aln->query->infer_CDR3_aa,
                    $aln->query->infer_CDR3_aa_length
                );

                foreach my $r (@regions) {
                    my $pred = "has_$r";
                    if ( $sub_r->$pred ) {
                        push @region_col,
                          (
                            $sub_r->$r,
                            $feature->alignments->$r->gaps,
                            $feature->alignments->$r->mismatches,
                            $sub_t->$r,
                          );
                    }
                    else {
                        push @region_col, ('N/A') x 4;
                    }
                }

                push @col, @region_col;

            }
            else {
                #say "No alignment to " . $feature->query_id;
                push @col, ('N/A') x ( scalar(@col_header) - 1 );
                $no_aln++;
            }
            $i++;

            say join "\t", @col;
        }
        $self->log->warn("==> END $cmd <==");
    }
}

class Main {
    import MyApp;
    MyApp->new_with_command->run();
}

=head1 NAME 

    MyApp

=head1 SYNOPSIS
  This application requires Perl 5.10.0 or higher   
  This application requires, at least, the following modules to work:
    - Moose
    - MooseX::App::Command

  Here, you want to concisely show a couple of SIMPLE use cases.  You should describe what you are doing and then write code that will run if pasted into a script.  

  For example:

  USE CASE: PRINT A LIST OF PRIMARY IDS OF RELATED FEATURES

    my $gene = new Modware::Gene( -feature_no => 4161 );

    foreach $feature ( @{ $gene->features() } ) {
       print $feature->primery_id()."\n";
    }

=head1 DESCRIPTION

   Here, AT A MINIMUM, you explain why the object exists and where it might be used.  Ideally you would be very detailed here. There is no limit on what you can write here.  Obviously, lesser used 'utility' objects will not be heavily documented.

   For example: 

   This object attempts to group together all information about a gene
   Most of this information is returned as references to arrays of other objects.  For example
   the features array is one such association.  You would use this whenever you want to read or write any 
   properties of a gene.


=head1 AUTHOR

Thiago Yukio Kikuchi Oliveira E<lt>stratust at gmail.comE<gt>

Copyright (c) 2013 Rockefeller University - Nussenzweig's Lab

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html

=head1 METHODS

=cut
