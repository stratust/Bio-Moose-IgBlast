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

class MyApp::Foo {
    extends 'MyApp'; # inherit log
    use MooseX::App::Command;    # important
    use MooseX::FileAttribute;
    use Bio::Moose::IgBlastIO;
    use Data::Printer deparse =>1, class => {expand => 'all'};
    use Data::Dumper;

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
        
        my $in = Bio::Moose::IgBlastIO->new( file => $self->input_file, format => 'format3' );

        while (my $feature = $in->next_feature){
            #say $feature->query_id;
            #p $feature->CDR3 if $feature->CDR3;
            #say $feature->domain_classification;
            #p $feature if $feature->query_id eq '3A_P95_G123_H1-IgGinternal.ab1';
            p $feature;
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

