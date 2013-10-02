use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::IgBlastIO {
    use Carp::Always;
    use Bio::Moose::IgBlastI::Format3;
    use Bio::Moose::IgBlastI::Format4;
    use Bio::Moose::IgBlastI::Format7;
    use MooseX::StrictConstructor;
  
    has 'file' => (
        is            => 'ro',
        isa           => 'Str | Path::Class::File',
        required      => 1,
        documentation => 'IgBlast file to be open',
    );
 
    has 'format' => (
        is            => 'ro',
        isa           => 'Str',
        required      => 1,
        documentation => 'Format of IgBlast to parse: format3, format4 or format7',
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
        my %index = (
            format3 => 'Bio::Moose::IgBlastI::Format3',
            format4 => 'Bio::Moose::IgBlastI::Format4',
            format7 => 'Bio::Moose::IgBlastI::Format7',
        );
        
        if ($index{$self->format}){
           my $o = $index{$self->format}->new(file => $self->file);
           return $o->features;
        }
        else{
            die "Unknown IgBlast format ". $self->format;
        }
    }
}
