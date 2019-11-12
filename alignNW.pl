#!/usr/bin/perl

use strict;
use warnings;
use autodie;

use Carp;

use Algorithm::NeedlemanWunsch;
use Bio::SeqIO;

my $sequence_file = shift or die; # assumed FASTA format

# Default primers are EMP 16S:
# http://www.earthmicrobiome.org/protocols-and-standards/16s/

my $fwd_primer = shift || q{GTGYCAGCMGCCGCGGTAA}; 
my $rev_primer = shift || q{GGACTACNVGGGTWTCTAAT}; 

die qq{no such file '$sequence_file'\n} if ! -f $sequence_file;

my $seqio = Bio::SeqIO->new( -file => $sequence_file, -format => 'fasta');

my $seqobj = $seqio->next_seq();

( my $dnaseq = $seqobj->seq() ) =~ s{ U }{T}gxms;

print qq{$dnaseq\n};

my $rev_primer_revcmp = reverse_complement($rev_primer);

print qq{forward primer:\t$fwd_primer\n},
      qq{reverse primer:\t$rev_primer\n},
      qq{reverse complement of reverse primer:\t$rev_primer_revcmp\n};

pairwise_align( sequence1   => $dnaseq,
                sequence2   => [ $fwd_primer, $rev_primer_revcmp ],
                local_align => 1,
                verbosity   => 0,
              );

#print qq{$score\n\n}, join(qq{\n}, @alignment), qq{\n};

exit 0;


sub create_subst_matrix {

    my $matrix_name = shift;

    my %matrix;

=pod
    The EMBOSS EDNAFULL matrix
    /usr/share/EMBOSS/data/EDNAFULL
    Note that the treatment of ambiguity codes is arguably unsatisfactory;
    it depends on whether you want to view the match as 'best case' or
    allow for less than best case.
=cut

    $matrix{ednafull} = <<END_EDNAFULL;
    A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N   U
A   5  -4  -4  -4  -4   1   1  -4  -4   1  -4  -1  -1  -1  -2  -4
T  -4   5  -4  -4  -4   1  -4   1   1  -4  -1  -4  -1  -1  -2   5
G  -4  -4   5  -4   1  -4   1  -4   1  -4  -1  -1  -4  -1  -2  -4
C  -4  -4  -4   5   1  -4  -4   1  -4   1  -1  -1  -1  -4  -2  -4
S  -4  -4   1   1  -1  -4  -2  -2  -2  -2  -1  -1  -3  -3  -1  -4
W   1   1  -4  -4  -4  -1  -2  -2  -2  -2  -3  -3  -1  -1  -1   1
R   1  -4   1  -4  -2  -2  -1  -4  -2  -2  -3  -1  -3  -1  -1  -4
Y  -4   1  -4   1  -2  -2  -4  -1  -2  -2  -1  -3  -1  -3  -1   1
K  -4   1   1  -4  -2  -2  -2  -2  -1  -4  -1  -3  -3  -1  -1   1
M   1  -4  -4   1  -2  -2  -2  -2  -4  -1  -3  -1  -1  -3  -1  -4
B  -4  -1  -1  -1  -1  -3  -3  -1  -1  -3  -1  -2  -2  -2  -1  -1
V  -1  -4  -1  -1  -1  -3  -1  -3  -3  -1  -2  -1  -2  -2  -1  -4
H  -1  -1  -4  -1  -3  -1  -3  -1  -3  -1  -2  -2  -1  -2  -1  -1
D  -1  -1  -1  -4  -3  -1  -1  -3  -1  -3  -2  -2  -2  -1  -1  -1
N  -2  -2  -2  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -2
U  -4   5  -4  -4  -4   1  -4   1   1  -4  -1  -4  -1  -1  -2   5
END_EDNAFULL

=pod

    This modified version takes the 'assume the best (usually)' approach.
    This makes sense if one of the sequences is a primer sequence (so in that
    context it represents more than one actual sequence).

    So not only does, e.g. W get a maximum score with both A and T, it also
    gets a maximum score with W. That would not be the correct score from a
    probabilistic sense, because there is effectively a 1/2 chance that the
    two bases are the same (assuming a completely uniform background frequency
    of all 4 bases, which is a caveat).

    However with overlapping ambiguity codes such
    as Y (C or T), W (A or T) has the same score (-2) as in the original.
    All such overlaps retain their original scores. Examples:

    S (G or C) versus B (C, G or T) has a score of -1; there is effectively a
    1/3 chance that the bases are literally the same.
    S (G or C) versus M (A or C) has a score of -2; there is effectively a
    1/4 chance that the bases are literally the same.
    S (G or C) versus H (A, C or T) has a score of -3; there is effectively a
    1/6 chance that the bases are literally the same.
    D (A, G or T) versus V (A, C or G) has a score of -2; there is effectively a
    2/9 chance that the bases are literally the same.

    If the matrix were really catering for the "one sequence is a primer"
    scenario, then again the above scores are not correct, but you would need
    to know which sequence is the primer. But that scenario would really mean
    that there could be an ambiguity code in only one sequence anyway, so it is
    moot. So I have not changed those ambiguity versus ambiguity scores.

    All scores involving N have been changed to the maximum. Again, that is
    not a correct probability in one sense (e.g. A v N has a 1/4 chance of
    being a real match); but in the primer sense, there is a 100% chance of a
    match.

    An alternative approach might be to have a maximum score only when one of
    the two bases compared is unambiguous. In that case, the score of W v W
    would remain as in the original (-1, instead of 5 here).

    This matrix below is intended to be used for comparing primer sequences
    (which often include ambiguity codes) versus known genomic DNA sequences
    of potential amplicons (which don't).

=cut

    $matrix{ednafull_liberal} = <<END_EDNAFULLLIB;
    A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N   U
A   5  -4  -4  -4  -4   5   5  -4  -4   5  -4   5   5   5   5  -4
T  -4   5  -4  -4  -4   5  -4   5   5  -4   5  -4   5   5   5   5
G  -4  -4   5  -4   5  -4   5  -4   5  -4   5   5  -4   5   5  -4
C  -4  -4  -4   5   5  -4  -4   5  -4   5   5   5   5  -4   5  -4
S  -4  -4   1   1   5  -4  -2  -2  -2  -2  -1  -1  -3  -3   5  -4
W   1   1  -4  -4  -4   5  -2  -2  -2  -2  -3  -3  -1  -1   5   1
R   1  -4   1  -4  -2  -2   5  -4  -2  -2  -3  -1  -3  -1   5  -4
Y  -4   1  -4   1  -2  -2  -4   5  -2  -2  -1  -3  -1  -3   5   1
K  -4   1   1  -4  -2  -2  -2  -2   5  -4  -1  -3  -3  -1   5   1
M   1  -4  -4   1  -2  -2  -2  -2  -4   5  -3  -1  -1  -3   5  -4
B  -4  -1  -1  -1  -1  -3  -3  -1  -1  -3   5  -2  -2  -2   5  -1
V  -1  -4  -1  -1  -1  -3  -1  -3  -3  -1  -2   5  -2  -2   5  -4
H  -1  -1  -4  -1  -3  -1  -3  -1  -3  -1  -2  -2   5  -2   5  -1
D  -1  -1  -1  -4  -3  -1  -1  -3  -1  -3  -2  -2  -2   5   5  -1
N   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5
U  -4   5  -4  -4  -4   1  -4   1   1  -4  -1  -4  -1  -1   5   5
END_EDNAFULLLIB

    die qq{unknown matrix '$matrix_name'\n} if !exists $matrix{lc $matrix_name};

    matrix_from_string($matrix{lc $matrix_name});

}

sub matrix_from_string {

    

}

sub reverse_complement {

    my $sequence = shift
      or croak qq{specify sequence for reverse complementing};

    my %rev_cmp_of = (

        A => 'T',
        G => 'C',
        C => 'G',
        T => 'A',
        Y => 'R',
        R => 'Y',
        W => 'S',
        S => 'W',
        K => 'M',
        M => 'K',
        D => 'H',
        V => 'B',
        H => 'D',
        B => 'V',
        U => 'A',
        N => 'N',
        X => 'N',

    );

    my $seqlen = length $sequence;

    my $revcmp = q{};

    for my $i ( reverse 0 .. $seqlen - 1 ) {
        my $base = uc substr $sequence, $i, 1;
        $revcmp .= $rev_cmp_of{$base};
    }

    return $revcmp;
}

sub pairwise_align {

=pod

    Uses the Algorithm::NeedlemanWunsch module to perform 1 or more pairwise
    sequence alignments.

    The first sequence is involved in all of the alignment(s), and is aligned
    in turn with all the others. If the client code needs to align a sequence
    with more than one other sequence, then it is more efficient to make a
    single call to pairwise_align(), i.e. of the form
    pairwise_align(sequence1 => $first_seq,
                   sequence2 => [ $second_seq, $third_seq,..])
    especially if the first sequence is long. One reason for that is that the
    input to Algorithm::NeedlemanWunsch->align() includes two arrayrefs,
    where each ref'd array contains all of the sequence letters of bases or
    residues as a list, i.e. one letter per element. With a single call of
    pairwise_align(), the first sequence need be prepared only once.

=cut

    my %arg = @_;

    my $sequence1 = delete $arg{sequence1}; # expects scalar
    my $sequence2 = delete $arg{sequence2}; # scalar or ref to array of scalars
    my $type      = delete $arg{type}; # 'dna', 'rna', 'protein'

    croak qq{no first sequence specified}    if !defined $sequence1;
    croak qq{no other sequence(s) specified} if !defined $sequence2;

    $type ||= 'dna';

    my $ref = ref $sequence2;

    my @sequences;

    if ($ref eq 'ARRAY') {
        @sequences = @{$sequence2};
    }
    else {
        @sequences = ( $sequence2 );
    }

    $sequence1 = uc $sequence1;
    $sequence1 =~ s{ U }{T}gxms if lc($type) eq 'dna';

    my @seqletters1 = split m{}xms, $sequence1;
#    my @seqletters_list;

    my $n_alignment = 0;
    for my $seq (@sequences) {
        $seq = uc $seq;
        $seq =~ s{ U }{T}gxms if lc($type) eq 'dna';
        my @seqletters = split m{}xms, $seq;
#        push @seqletters_list, \@seqletters;
        my @alignment = 
            NWalign(seqlist1 => \@seqletters1,
                seqlist2 => \@seqletters ,
                %arg);
        my $score = pop @alignment;
        print qq{alignment }, ++$n_alignment, qq{:\n\n},
            qq{$alignment[0]\n$alignment[1]\n\nscore: $score\n\n};
    }

}

sub NWalign {

    my %arg = @_;

    my ($seq_aref1, $seq_aref2, $gop, $gep, $local_align, $verbosity) =
        @arg{'seqlist1', 'seqlist2', 'gap_open',
             'gap_extend', 'local_align', 'verbosity'};

    my @seq_aref = ($seq_aref1, $seq_aref2);

    $verbosity   ||= 0;
    $local_align ||= 0;

    my $matcher = Algorithm::NeedlemanWunsch->new(\&score_sub);

    if ( defined($gop) || defined($gep) ) {
        croak qq{both or neither GOP and GEP must be defined\n}
          if !(defined($gop) && defined($gep));
        # Very surprisingly, the Algorithm::NeedlemanWunsch->align() function
        # requires that the gap opening penalty is **less** than the gap
        # extension penalty. Nonetheless, the default values seem to work well
        # with the tests I've done of long sequences versus very short (primer)
        # sequences.
        $matcher->gap_open_penalty($gop);
        $matcher->gap_extend_penalty($gep);
    }

    $matcher->local($local_align); # boolean

    my @alignment = ( q{}, q{} );

    # Because callbacks provide no means of passing additional arguments
    # nor of returning values, coderefs are defined here, which access the
    # variables in the current subroutine's scope.

    my $on_align_cref = sub {
        my @arg = @_; # there are two values passed to this callback:
                      # the indices of the aligned letter in the first and
                      # second sequences, respectively
        my $i = 0;

        for my $arg (@arg) {
            print qq{on_align\t}, $i, qq{, '$arg', },
                  $seq_aref[$i]->[$arg], qq{\n} if $verbosity;
            $alignment[$i] = $seq_aref[$i]->[$arg]. $alignment[$i];
            $i++;
        }
        print qq{\n} if $verbosity;
    };

    my $on_shift_a_cref = sub {
        my $i = 0;
        for my $arg (@_) {
            print qq{on_shift_a\t}, $i, qq{, $arg, },
                  $seq_aref[0]->[$arg], qq{\n\n} if $verbosity;
            $alignment[0] = $seq_aref[0]->[$arg]. $alignment[0];
            $alignment[1] = q{-} . $alignment[1];
            $i++;
        }
    };

    my $on_shift_b_cref = sub {
        my $i = 0;
        for my $arg (@_) {
            print qq{on_shift_b\t}, $i, qq{, $arg, },
                  $seq_aref[1]->[$arg], qq{\n\n} if $verbosity;
            $alignment[0] = q{-} . $alignment[0];
            $alignment[1] = $seq_aref[1]->[$arg]. $alignment[1];
            $i++;
        }
    };

    my $score = $matcher->align(
        $seq_aref[0],
        $seq_aref[1],
        {
            align   => $on_align_cref  ,
            shift_a => $on_shift_a_cref,
            shift_b => $on_shift_b_cref,
        }
    );

    return @alignment, $score;
}

sub score_sub {
    if (!@_) {
        return -2; # gap penalty
    }
 
    return ($_[0] eq $_[1]) ? 5 : -4;
}

=pod 
sub on_align {

    my @arg = @_;
    my $i = 0;

    for my $arg (@arg) {
        print qq{on_align\t}, $i, qq{, '$arg', }, $seq[$i]->[$arg], qq{\n};
        $alignment[$i] = $seq[$i]->[$arg]. $alignment[$i];
        $i++;
    }
    print qq{\n};
}

sub on_shift_a {

    my $i = 0;
    for my $arg (@_) {
        print qq{on_shift_a\t}, $i, qq{, $arg, }, $seq[0]->[$arg], qq{\n};
        $alignment[0] = $seq[0]->[$arg]. $alignment[0];
        $alignment[1] = q{-} . $alignment[1];
        $i++;
    }

}

sub on_shift_b {

    my $i = 0;
    for my $arg (@_) {
        print qq{on_shift_b\t}, $i, qq{, $arg, }, $seq[1]->[$arg], qq{\n};
        $alignment[0] = q{-} . $alignment[0];
        $alignment[1] = $seq[1]->[$arg]. $alignment[1];
        $i++;
    }

}

=cut


#my @aligned = ( [ ], [ ] );


