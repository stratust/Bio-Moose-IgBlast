use MooseX::Declare;
use Method::Signatures::Modifiers;

role MyType {
    use Moose::Util::TypeConstraints;
    use Bio::Seq;

    subtype 'BioSeq' => as "Bio::Seq";

    coerce 'BioSeq', 
        from 'Str', 
        via {
            my $value = shift;
            my $object = Bio::Seq->new(
                    -seq      => $value,
                    -alphabet => 'dna'
                    );
            return $object;
        };
}

class Bio::Moose::IgBlast::HitTable::Entry {
    with 'MyType';
    use MooseX::StrictConstructor;
    use MooseX::Aliases;
    
    has 'query_id'         => ( is => 'ro', isa => 'Str', required => 1, alias => 'qseqid' );
    has 'subject_id'       => ( is => 'ro', isa => 'Str', required => 1, alias => 'sseqid' );
    has 'identity'         => ( is => 'ro', isa => 'Num', required => 1, alias => 'pident' );
    has 'alignment_length' => ( is => 'ro', isa => 'Int', required => 1, alias => 'length' );
    has 'mismatches'       => ( is => 'ro', isa => 'Int', required => 1, alias => 'mismatch' );
    has 'gaps'             => ( is => 'ro', isa => 'Int', required => 0 );
    has 'gap_opens'        => ( is => 'ro', isa => 'Int', required => 1, alias => 'gapopen' );
    has 'query_start'      => ( is => 'ro', isa => 'Int', required => 1, alias => 'qstart' );
    has 'query_end'        => ( is => 'ro', isa => 'Int', required => 1, alias => 'qend' );
    has 'subject_start'    => ( is => 'ro', isa => 'Int', required => 1, alias => 'sstart' );
    has 'subject_end'      => ( is => 'ro', isa => 'Int', required => 1, alias => 'send' );
    has 'evalue'           => ( is => 'ro', isa => 'Num', required => 1 );
    has 'bit_score'        => ( is => 'ro', isa => 'Num', required => 1, alias => 'bitscore' );

    has [qw/sseq qseq/] => ( is => 'ro', isa => 'BioSeq', coerce => 1 );

    
    has [
        qw/btop frames nident positive ppos qacc qaccver qframe qgi qlen sacc saccver sallacc sallgi sallseqid score sframe sgi slen /
    ] => ( is => 'ro', isa => 'Str' );
}
