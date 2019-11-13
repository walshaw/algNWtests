#!/usr/bin/perl

# 
use strict;
use warnings;
use autodie;

use Carp;

use Getopt::Long;

use Data::Dumper;
use List::Util      qw( first min max );
use List::MoreUtils qw( first_index uniq );

use Algorithm::NeedlemanWunsch;
use Bio::SeqIO;

my $sequence_file;
my $sequence_format = q{fasta};
my $fwd_primer = q{GTGYCAGCMGCCGCGGTAA};  # EMP 16S 515F
my $rev_primer = q{GGACTACNVGGGTWTCTAAT}; # EMP 16S 806R
my $matrix_name;
my $local_alignment;  # boolean; false = global, true = local; this refers
                      #     to the variant of the algorithm applied by
                      #     Algorithm::NeedlemanWunsch; *not* to how the
                      #     alignment is fetched/shown (if at all);
                      #     even if the alignment algorithm is local, you
                      #     might still want to view the location of (say)
                      #     a short sequence aligned to the whole of a long
                      #     sequence
###my $create_alignment; # boolean; false = no alignment created or shown
my $create_alignment; # 0 or undef means no alignment; other options are
                      # 1 trim (if necessary) alignment so that it spans
                      # sequence 1 and no further;
                      # 2 trim (if necessary) alignment so that it spans
                      # sequence 2 and no further;
                      # 3 trim such that the alignment begins with the
                      # first no-gap position and ends with the last
                      # no-gap position (this is the harshest trim,
                      # effectively producing a traditional local-alignment);
                      # 4 no trimming;
                      # instead of the numbers, can use 'none' for 0;
                      # 'trimmed1' for 1; 'trimmed2' for 2; 'trimmed' for 3;
                      # 'full', 'all', 'global' for 4.
my $get_coords1;      # boolean; refer to arg name 'coords1' in NWalign()
my $get_coords2;      # boolean; refer to arg name 'coords2' in NWalign()
my $verbosity = 0;

my @alignment_type = qw( none trimmed1 trimmed2 trimmed global );

my $usage = qq{Usage:\n\n$0 [ -sequences ] SEQUENCE_FILE [ options ]\n};

die $usage if !GetOptions(

    "sequences|sequence_file|file=s" => \$sequence_file,
    "format=s"                       => \$sequence_format,
    "forward_primer|fwd_primer=s"    => \$fwd_primer,
    "reverse_primer|rev_primer=s"    => \$rev_primer,
    "matrix=s"                       => \$matrix_name,
    "local_scoring|local_alignment"  => \$local_alignment, # better as 'local_scoring'?
    "alignment=s"                    => \$create_alignment,# then this can have value 'local'?
    "coords1"                        => \$get_coords1,
    "coords2"                        => \$get_coords2,
    "verbosity=i"                    => \$verbosity,
);

$sequence_file ||= shift or die $usage;

$create_alignment ||= 0;

if ($create_alignment =~ m{ \A \d \z }xms) {
    die qq{invalid alignment type (should be 0 .. 4)\n}
      if ($create_alignment < 0) || ($create_alignment > 4);
}
else {
    $create_alignment =~ s{ \A (?:full|all) .* \z }{global}ixms;
    my $idx = first_index { lc $_ eq $create_alignment } @alignment_type;
    die qq{invalid alignment type '$create_alignment'} if $idx == -1;
    $create_alignment = $idx;
}

print qq{alignment type: $alignment_type[$create_alignment]\n};


# Default primers are EMP 16S:
# http://www.earthmicrobiome.org/protocols-and-standards/16s/

###my $matrix_name = shift; #  || q{ednafull};

die qq{no such file '$sequence_file'\n} if ! -f $sequence_file;



my $seqio = Bio::SeqIO->new( -file   => $sequence_file,
                             -format => $sequence_format);

my $seqobj = $seqio->next_seq();

( my $dnaseq = $seqobj->seq() ) =~ s{ U }{T}gxms;

print qq{$dnaseq\n};

my $rev_primer_revcmp = reverse_complement($rev_primer);

print qq{forward primer:\t$fwd_primer\n},
      qq{reverse primer:\t$rev_primer\n},
      qq{reverse complement of reverse primer:\t$rev_primer_revcmp\n};

# To ensure that the substitution matrix is created once and only once, it is
# done now, and then passed to the aligner

# will be undef if there is no $matrix_name
my $matrix_href = create_subst_matrix($matrix_name);

pairwise_align( sequence1   => $dnaseq,
                sequence2   => [ $fwd_primer, $rev_primer_revcmp ],
                local_align => $local_alignment,
                matrix      => $matrix_href, # can be undef
                alignment   => $create_alignment, 
                coords1     => $get_coords1, # can be undef
                coords2     => $get_coords2, # can be undef
                verbosity   => 0,
              );

#print qq{$score\n\n}, join(qq{\n}, @alignment), qq{\n};

exit 0;


sub create_subst_matrix {

    my $matrix_name = shift;

    return undef if !$matrix_name;

    my $matrices_href = set_matrix_strings();

    die qq{unknown matrix '$matrix_name'\n}
      if !exists $matrices_href->{lc $matrix_name};

    my $matrix_aref = 
        matrix_from_string($matrices_href->{lc $matrix_name});

    return $matrix_aref;
}

sub set_matrix_strings {

    # creates contents of this hash; values are long strings specifying a
    # substitution matrix

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
    as Y (C or T), W (A or T) the same score (-2) as in the original is us.
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
S  -4  -4   5   5   5  -4  -2  -2  -2  -2  -1  -1  -3  -3   5  -4
W   5   5  -4  -4  -4   5  -2  -2  -2  -2  -3  -3  -1  -1   5   5
R   5  -4   5  -4  -2  -2   5  -4  -2  -2  -3  -1  -3  -1   5  -4
Y  -4   5  -4   5  -2  -2  -4   5  -2  -2  -1  -3  -1  -3   5   5
K  -4   5   5  -4  -2  -2  -2  -2   5  -4  -1  -3  -3  -1   5   5
M   5  -4  -4   5  -2  -2  -2  -2  -4   5  -3  -1  -1  -3   5  -4
B  -4   5   5   5  -1  -3  -3  -1  -1  -3   5  -2  -2  -2   5   5
V   5  -4   5   5  -1  -3  -1  -3  -3  -1  -2   5  -2  -2   5  -4
H   5   5  -4   5  -3  -1  -3  -1  -3  -1  -2  -2   5  -2   5   5
D   5   5   5  -4  -3  -1  -1  -3  -1  -3  -2  -2  -2   5   5   5
N   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5   5
U  -4   5  -4  -4  -4   5  -4   5   5  -4   5  -4   5   5   5   5
END_EDNAFULLLIB

    return \%matrix;

}

sub matrix_from_string {

    my $matrix_string = shift or croak qq{no string supplied};

    my @rows = split m{ \n }xms, $matrix_string;

    print scalar @rows, qq{ rows in matrix (including header)\n};

    my @column_headings;
    my @row_headings;
    my %score; # the resulting matrix

    for my $row_string (@rows) {

        next if $row_string !~ m{ \S }xms;

        $row_string =~ s{ ( \A \s+ | \s+ \z ) }{}gxms;

        my @fields = split m{ \s+ }xms, $row_string;

        if (!@column_headings) {
            @column_headings = map { uc $_ } @fields;
            next;
        }

        croak sprintf qq{Expected %d columns including row header, }.
                      qq{but there are %d},
                scalar(@column_headings) + 1, scalar @fields
            if scalar(@fields) != scalar(@column_headings)+1;

        my $row_heading = uc (shift @fields);
        push @row_headings, $row_heading;

        for my $i ( 0 .. $#fields ) {

            # check for symmetry

            croak qq{non-symmetrical scores in matrix:\n}.
                  qq{    $column_headings[$i] v $row_heading = }.
                  $score{$column_headings[$i]}{$row_heading}.
                  qq{\n    $row_heading v $column_headings[$i] = $fields[$i]\n}
              if (defined($score{$column_headings[$i]}{$row_heading}) &&
                 ($score{$column_headings[$i]}{$row_heading} != $fields[$i]));

            $score{$row_heading}{$column_headings[$i]} = $fields[$i];

        }

    }

    # sanity-check row and column headings

    my @sorted_col_heading = sort @column_headings;
    my @sorted_row_heading = sort @row_headings;

    my @uniq_col_heading = uniq @sorted_col_heading;
    my @uniq_row_heading = uniq @sorted_row_heading;

   croak qq{non-unique letter keys in matrix column headings:},
        join(q{, }, @uniq_col_heading), qq{\n}
       if @uniq_col_heading != @sorted_col_heading;

   croak qq{non-unique letter keys in matrix row headings:},
        join(q{, }, @uniq_row_heading), qq{\n}
       if @uniq_row_heading != @sorted_row_heading;


    for my $i ( 0 .. $#sorted_col_heading ) {
        croak qq{column heading '$sorted_col_heading[$i]' }.
              qq{doesn't match row heading '$sorted_row_heading[$i]'}
          if $sorted_col_heading[$i] ne $sorted_row_heading[$i];
    }

 
    print qq{successfully read scoring matrix for },
      join(q{, }, @column_headings), qq{\n};

    return \%score;

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

sub match_count {

    # The two sequence strings must already be aligned.
    # They should be the same length (might include gap characters).
    # But if one is longer than the other, gap characters are assumed.
    my $aligned_seq1 = shift;
    my $aligned_seq2 = shift;
    my $optimistic   = shift;
    my $gap_char     = shift || q{-};

    # unoptimistic not yet supported (would calculate probabilities)
    croak qq{only 'optimistic' match counts supported}
      if !$optimistic;

    my %bases = (
        A => [ 'A' ],
        G => [ 'G' ],
        C => [ 'C' ],
        T => [ 'T' ],
        Y => [ 'C', 'T' ],
        R => [ 'A', 'G' ],
        W => [ 'A', 'T' ],
        S => [ 'G', 'C' ],
        K => [ 'G', 'T' ],
        M => [ 'A', 'C' ],
        D => [ 'A', 'G', 'T' ],
        V => [ 'A', 'C', 'G' ],
        H => [ 'A', 'C', 'T' ],
        B => [ 'C', 'G', 'T' ],
        U => [ 'T' ],
        N => [ 'A', 'C', 'G', 'T' ],
        X => [ 'A', 'C', 'G', 'T' ],
        $gap_char => [], # won't match anything, not even another gap
    );

    for my $b (keys %bases) {
        push @{$bases{$b}}, 'U' if first { $_ eq 'T' } @{$bases{$b}};
    }

    my $aln_length = max length($aligned_seq1), length($aligned_seq2);
    # show a warning if not the same length...

    for my $seq ($aligned_seq1, $aligned_seq2) {
        my $short = $aln_length - length $seq;
        $seq .= $gap_char x $short;
    }

    my $n_matches = 0;

    for my $i ( 0 .. $aln_length - 1 ) {

        my $base1 = substr $aligned_seq1, $i, 1;
        my $base2 = substr $aligned_seq2, $i, 1;

        my $n_options1 = scalar @{$bases{$base1}};
        my $n_options2 = scalar @{$bases{$base2}};
        my @union = uniq (@{$bases{$base1}}, @{$bases{$base2}});
        # If the two sets are mutually exclusive, then the size of the union
        # is the sum of the sizes of the two sets; otherwise, they have at
        # least one base in common, which is a 'match' in optimistic-speak
        $n_matches++ if scalar(@union) < $n_options1 + $n_options2;
    }

    return $n_matches, $n_matches / $aln_length;
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

    Some arguments (values of %arg, see below) are consumed by this current
    subroutine, i.e. with keys:

            'sequence1', 'sequence2' (see above)
            'dna'

    Most arguments are passed straight to NWalign, without any intervention
    in the current subroutine, by passing what remains of %arg (see below);
    keys are:

            'gap_open'   , 'gap_extend',
            'local_align',
            'matrix'     ,
            'alignment'  ,
            'coords1'    , 'coords2'   ,
            'verbosity'

    'alignment' is an exception; it is passed as a boolean to NWalign (directs
    that to either not create, or create, an alignment) but within the
    current subroutine pairwise_align() the different true values ( 1 .. 4)
    have different meanings, which are applied here (how to trim the alignment,
    if at all).

    The current subroutine also creates addition key => value pairs, with keys:

            'seqlist1'   , 'seqlist2'  ,


=cut

    my %arg = @_;

    my $sequence1 = delete $arg{sequence1}; # expects scalar
    my $sequence2 = delete $arg{sequence2}; # scalar or ref to array of scalars
    my $type      = delete $arg{type}; # 'dna', 'rna', 'protein'
    my $alignment_type = $arg{alignment}; # intentionally not deleted

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

        my ($alignment_score, %results) = 
            NWalign(
                seqlist1 => \@seqletters1,
                seqlist2 => \@seqletters ,
                %arg);
###        my $score = pop @alignment;
###        print qq{alignment }, ++$n_alignment, qq{:\n\n},

        print qq{score: $alignment_score:\n\n};

        if ($alignment_type) {

            my $display_alignment_aref =
              trim_alignment($results{alignment},$alignment_type);

            print qq{$display_alignment_aref->[0]\n},
                  qq{$display_alignment_aref->[1]\n\n};
            my ($matches, $matches_prpn) = match_count(@{$display_alignment_aref},1);
print qq{match count: $matches (}, 100 * $matches_prpn  , qq{%)\n};
            #print qq{$results{alignment}->[0]\n$results{alignment}->[1]\n\n}
            #  if $results{alignment};
        }

        show_coords($results{coords1}, scalar(@seqletters1),
            qq{coords of sequence 2 corresponding to coords of sequence 1},
            qq{seq1}, qq{seq2})
          if $results{coords1};

        show_coords($results{coords2}, scalar(@seqletters),
            q{coords of sequence 1 corresponding to coords of sequence 2},
            q{seq2}, q{seq1})
          if $results{coords2};

        # Save the results (conditionally requested) here, in a growing
        # array (one element per alignment); each element will include,
        # conditionally, the full (and display) alignment; the
        # table(s) of counterpart-corresponding coordinates; number of
        # mismatches in the display alignment, and thus the corresponding
        # pairwise identity
        # The alignment score will always be a result.

    }

}

sub show_coords {

    my $coords_aref = shift;
    my $n_coords    = shift || scalar(@{$coords_aref});
    my $title = shift || q{};
    my $heading1 = shift || q{1};
    my $heading2 = shift || q{1};

    printf qq{$title\n\n\t%8s\t%8s\n\n}, $heading1, $heading2;

    for my $i ( 0 .. $n_coords - 1 ) {
        my $paired_coord = $coords_aref->[$i];
        $paired_coord = q{-} if !defined $paired_coord;
        printf qq{\t%8d\t%8s\n}, $i, $paired_coord;
    }

}

sub trim_alignment {

    my $alignment_aref = shift; # 2 elements; 1 string per sequence
    my $alignment_type = shift || 4; # default just returns the original

=pod

    The meaning of $alignment_type:

    Notionally, 0 means that no alignment results; but if that's the case,
    then this subroutine is unlikely to have been called in the first place.

    1: The alignment is trimmed at both ends, such that the resulting
    alignment begins with the first letter of the first sequence, and ends
    with the last letter of the first sequence.

    2: The alignment is trimmed at both ends, such that the resulting
    alignment begins with the first letter of the second sequence, and ends
    with the last letter of the second sequence.

    3: The alignment is trimmed at both ends, such that the resulting
    alignment begins with the first position that is occupied by a letter
    (rather than a gap character) of both sequences, and ends with the last
    position that is occupied by a letter of both sequences.

    4: The alignment is not trimmed; in other words, the whole of both
    sequences is included; note that this is the case even if the local-
    alignment variant of the algorith is applied (see -local_alignment).
    The local alignment algorithm simply affects the way that gaps at
    each end of the aligned sequences are scored.

    For example, if you always want the returned alignment to represent a
    local alignment in the traditional sense - i.e. parts of both sequences
    which do not align with the other sequence are trimmed off (but only at
    both ends), then that is type 3. E.g.:

    sequence 1  ACGCGCTTTC--CGCGAG-----
    sequence 2  -----CTTTCGAGGAGAGTAATA

    becomes with type 3:

    sequence 1  CTTTC--CGCGAG
    sequence 2  CTTTCGAGGAGAG

    whereas with type 1:

    sequence 1  ACGCGCTTTC--CGCGAG
    sequence 2  -----CTTTCGAGGAGAG

    with type 2:

    sequence 1  CTTTC--CGCGAG-----
    sequence 2  CTTTCGAGGAGAGTAATA

    If you are expecting one sequence to be much shorter than the other, and
    you want a representation of to where on the longer sequence the
    shorter one aligns; e.g. if the first sequence was a long genomic sequence
    and the second represented a short PCR primer; then use (in this example)
    type 2.

    Note that the trimming positions are not determined from the tables of
    counterpart-aligned coordinates of each sequence - because the creation of
    those tables is only conditional.

=cut

    return if !$alignment_type; # maybe returning [ q{}, q{} ] is better?

    return $alignment_aref if $alignment_type == 4;

    my $aln_first;
    my $aln_last;

    my $first1;
    my $last1;
    my $first2;
    my $last2;

#    my ($first1, $last1) = non_gap_limits($alignment_aref->[0]);
#    my ($first2, $last2) = non_gap_limits($alignment_aref->[1]);

    ALN_TYPE: {

        (($alignment_type == 1) || ($alignment_type == 3)) && do {

            ($first1, $last1) = non_gap_limits($alignment_aref->[0]);

            ($alignment_type == 1) && do {
                $aln_first = $first1;
                $aln_last  = $last1;
                last ALN_TYPE;
            };

        };

        (($alignment_type == 2) || ($alignment_type == 3)) && do {
            ($first2, $last2) = non_gap_limits($alignment_aref->[1]);

            ($alignment_type == 2) && do {
                $aln_first = $first2;
                $aln_last  = $last2;
                last ALN_TYPE;
            };
        };

        ($alignment_type == 3) && do {

            $aln_first = max $first1, $first2;
            $aln_last  = min $last1 , $last2;
            last ALN_TYPE;
        };

        croak qq{invalid alignment-trimming code: $alignment_type};

    } # END block ALN_TYPE

    my $trimmed_length = 1 + $aln_last - $aln_first;

    my @trimmed_alignment =
      map { substr $_, $aln_first, $trimmed_length } @{$alignment_aref};

    return \@trimmed_alignment;

}

sub non_gap_limits {

    my $string   = shift;
    my $gap_char = shift || q{-};

    my $length = length $string;

    my $first_non_gap;
    my $last_non_gap;

    # Maybe more efficient to do this with a s/// substitution and then use
    # the lengths of the bits trimmed off to determine the limit?

    for my $i ( 0 .. $length - 1) {

        if (!defined($first_non_gap) && substr($string, $i, 1) ne $gap_char) {
            $first_non_gap = $i;
            return ($first_non_gap, $last_non_gap) if defined $last_non_gap;
        }

        if (!defined($last_non_gap) &&
            substr($string, -1*($i+1), 1) ne $gap_char) {
            $last_non_gap = $length - $i - 1;
            return ($first_non_gap, $last_non_gap) if defined $first_non_gap;
        }

    }

    return ($first_non_gap, $last_non_gap); # can only be (undef, undef)

}

sub NWalign {

    my %arg = @_;

    my ($seq_aref1     , $seq_aref2,
        $gop           , $gep,
        $local_align   , # boolean; if true, a local alignment is performed
        $matrix_href   ,
###        $alignment_aref, # optional output; ref to array of 2 string elements
###        $coords1_aref  , $coords2_aref, # both optional output; see below
        $create_alignment, # boolean
        $compile_coords1  , $compile_coords2,
        $verbosity) =

        @arg{'seqlist1'   , 'seqlist2'  ,
             'gap_open'   , 'gap_extend',
             'local_align',
             'matrix'     ,
             'alignment'  ,
             'coords1'    , 'coords2'   ,
             'verbosity'};

    my $no_actions = !($verbosity || $create_alignment ||
                       $compile_coords1 || $compile_coords2);

###    my @create_coords = ($create_coords1, $create_coords2);

=pod

    If the value of the coords1 key is true, then
    an array is created, with the same number of elements as there are letters
    in sequence 1; each value is the coordinate of sequence2 to which the
    letter in sequence1 has been aligned. To repeat: the indices represent the
    coords of sequence1 (corresponding to the indices of @{$seq_aref1}) and
    the values represent the coords of sequence2 (corresponding to the indices
    of @{$seq_aref2}).

    A ref to that array is then returned as the value of the 'coords1' key of
    the %optional_results hash.

    Conversely, if coords2 has a true value an array
    is created, with the same number of elements as there are letters in
    sequence 2; each value is the coordinate of sequence1 to which the letter
    in sequence2 has been aligned. A ref to that array is then returned as the
    value of the 'coords2' key of the %optional_results hash.

    These arrays can be requested optionally and independently for the sake of
    efficiency; if NWalign is being called many times, then building these
    arrays could be costly if they are never used.

    Usage example: sequence1 is long, e.g. a gene sequence (such as the 16S
    rRNA gene); sequence2 is a short sequence, representing a primer. You
    might then want to create the coords2 array, which will, for each
    coordinate of this short primer, contain the corresponding coordinate in
    the gene sequence; in other words, where in the gene the amplicon would
    begin (including the primer-matching segment itself). The coords1 array
    would be less useful, and the vast majority of its elements would have
    an undefined value.

=cut

    my @seq_aref = ($seq_aref1, $seq_aref2);

    $verbosity   ||= 0;
    $local_align ||= 0;
    $gop         ||= -5;

    my $score_cref;

    # One of two alternative forms of scoring function coderef is
    # defined here; the functionality of the scoring function dictated by
    # the module is a little cryptic. If it's called with no arguments,
    # then it must return "the gap penalty"; if it's called with 2 arguments,
    # then these are the two letters aligned together at that point in the
    # alignment. In which case, naturally it returns the score of the
    # alignment of those two letter.
    # So when is it called with no arguments? Only by the
    # constructor, if the calling of the constructor itself has only 1
    # argument. If the constructor is called with 2 arguments, then the
    # second argument is "the gap penalty" - thus, there is no need to ever
    # call the scoring function with no arguments.
    # If you've followed all that, then the next question is, how does "the
    # gap penalty" relate to the affine penalties, GOP and GEP (defined
    # respectively by calling the options gap_open_penalty() and
    # gap_extend_penalty() functions)? I'm not absolutely sure - see notes
    # on this below, but the module doc does imply that the 2-penalty
    # approach is an alternative to the default 1-penalty approach.
    # So the logic applied here is:
    #   no penalties passed to this function (NWalign())
    #       - then set the GOP to a default (see above) and explicitly
    #         specify it as arg 2 when the constructor is called; thus
    #         there is no need to implement a zero-arg mode in the scoring
    #         function;
    #   GOP passed only
    #       - then treat this in the same way as the default GOP, above;
    #   GEP passed only
    #       - then set the GOP to the default and use both; because neither
    #         the gap_open_penalty() nor gap_extend_penalty() functions
    #         may be called unless *both* are (see below)
    #   GOP and GEP passed
    #       - then call both gap_open_penalty() and gap_extend_penalty()
    

    if (defined $matrix_href) {

       $score_cref = sub {

           croak qq{'$_[0]' not known in scoring matrix}
             if !defined $matrix_href->{$_[0]};

           croak qq{'$_[1]' not known in scoring matrix}
             if !defined $matrix_href->{$_[0]}->{$_[1]};

           return $matrix_href->{$_[0]}->{$_[1]};

       };

    }
    else {

        $score_cref = sub {
            # effectively the default scoring matrix
            return ($_[0] eq $_[1]) ? 5 : -4;
        };

    }

    my $matcher = Algorithm::NeedlemanWunsch->new($score_cref, $gop);

    if ( defined $gep ) {
#    if ( defined($gop) || defined($gep) ) {
#        croak qq{both or neither GOP and GEP must be defined\n}
#          if !(defined($gop) && defined($gep));
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
    my @coords1;
    my @coords2;

    # Because callbacks provide no means of passing additional arguments
    # nor of returning values, coderefs are defined here, which access the
    # variables in the current subroutine's scope.

    my $on_align_cref = sub {

        return if $no_actions;

        my @coords = @_; # there are two values passed to this callback:
                      # the indices of the aligned letter in the first and
                      # second sequences, respectively
        my $i = 0;

        for my $coord (@coords) {

            print qq{on_align\t}, $i, qq{, '$coord', },
                  $seq_aref[$i]->[$coord], qq{\n} if $verbosity;

            $alignment[$i] = $seq_aref[$i]->[$coord] . $alignment[$i]
              if $create_alignment;

            $i++;
        }

        $coords1[$coords[0]] = $coords[1] if $compile_coords1;
        $coords2[$coords[1]] = $coords[0] if $compile_coords2;

        print qq{\n} if $verbosity;
    };

    my $on_shift_a_cref = sub {

        return if $no_actions;

        my $i = 0;

        for my $arg (@_) {

            print qq{on_shift_a\t}, $i, qq{, $arg, },
                  $seq_aref[0]->[$arg], qq{\n\n} if $verbosity;

            if ($create_alignment) {
                $alignment[0] = $seq_aref[0]->[$arg] . $alignment[0];
                $alignment[1] = q{-} . $alignment[1];
            }

            $i++;
        }
    };

    my $on_shift_b_cref = sub {

        return if $no_actions;

        my $i = 0;

        for my $arg (@_) {

            print qq{on_shift_b\t}, $i, qq{, $arg, },
                  $seq_aref[1]->[$arg], qq{\n\n} if $verbosity;

            if ($create_alignment) {
                $alignment[0] = q{-} . $alignment[0];
                $alignment[1] = $seq_aref[1]->[$arg]. $alignment[1];
            }

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

    my %optional_results;
    $optional_results{alignment} = \@alignment if $create_alignment;
    $optional_results{coords1}   = \@coords1   if $compile_coords1;
    $optional_results{coords2}   = \@coords2   if $compile_coords2;

    return $score, %optional_results; 
}

=head1 AUTHOR

John Walshaw


=head1 LICENCE AND COPYRIGHT

Copyright (c) 2019, John Walshaw. All rights reserved.

Developed in my own time using my own resources.

This is free software; you can redistribute it and/or
modify it under the same terms as Perl itself. See L<perlartistic>.


=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.


