#!/usr/bin/env perl
#
#              INGLÊS/ENGLISH
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  http://www.gnu.org/copyleft/gpl.html
#
#
#             PORTUGUÊS/PORTUGUESE
#  Este programa é distribuído na expectativa de ser útil aos seus
#  usuários, porém NÃO TEM NENHUMA GARANTIA, EXPLÍCITAS OU IMPLÍCITAS,
#  COMERCIAIS OU DE ATENDIMENTO A UMA DETERMINADA FINALIDADE.  Consulte
#  a Licença Pública Geral GNU para maiores detalhes.
#  http://www.gnu.org/copyleft/gpl.html
#
#  Copyright (C) 2019  Universidade Estadual Paulista "Júlio de Mesquita Filho"
#
#  Universidade Estadual Paulista "Júlio de Mesquita Filho" (UNESP)
#  Faculdade de Ciências Agrárias e Veterinárias (FCAV)
#  Laboratório de Bioinformática (LB)
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
#  http://www.fcav.unesp.br
#
# $Id$

=head1 NAME

=head1 SYNOPSIS

=head1 ABSTRACT

=head1 DESCRIPTION
    
    Arguments:

        -h/--help   Help
        -l/--level  Log level [Default: FATAL] 
            OFF
            FATAL
            ERROR
            WARN
            INFO
            DEBUG
            TRACE
            ALL

=head1 AUTHOR

Daniel Guariz Pinheiro E<lt>dgpinheiro@gmail.comE<gt>

Copyright (C) 2019 Universidade Estadual Paulista "Júlio de Mesquita Filho"

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html


=cut

use strict;
use warnings;
use Readonly;
use Getopt::Long;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $left, $right, $single_left, $single_right);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "pe1|pe_left=s"=> \$left,
            "pe2|pe_right=s"=> \$right,
            "se1|se_left=s"=> \$single_left,
            "se2|se_right=s"=> \$single_right
    ) or &Usage();


if ($level) {
    my %LEVEL = (   
    'OFF'   =>$OFF,
    'FATAL' =>$FATAL,
    'ERROR' =>$ERROR,
    'WARN'  =>$WARN,
    'INFO'  =>$INFO,
    'DEBUG' =>$DEBUG,
    'TRACE' =>$TRACE,
    'ALL'   =>$ALL);
    $LOGGER->logdie("Wrong log level ($level). Choose one of: ".join(', ', keys %LEVEL)) unless (exists $LEVEL{$level});
    Log::Log4perl->easy_init($LEVEL{$level});
}


use HTML::Template;
use File::Basename;
use Cwd qw(abs_path realpath fast_abs_path);

use constant MAX_RD_LEN=>151;
use constant AVG_INS=>500;
use constant RD_LEN_CUTOFF=>100;

my @left_file=split(/,/, $left);
my @right_file=split(/,/, $right);
my @single_left_file=split(/,/, $single_left);
my @single_right_file=split(/,/, $single_right);

$LOGGER->logdie("The numbers of elements in left and right input files differ") if (scalar(@left_file)!=scalar(@right_file));

if (scalar(@single_left_file)>0) {
    $LOGGER->logdie("The numbers of elements in left and right paired-end files and sigle left input file differ") if (scalar(@left_file)!=scalar(@single_left_file));
}
if (scalar(@single_right_file)>0) {
    $LOGGER->logdie("The numbers of elements in left and right paired-end files and sigle right input file differ") if (scalar(@left_file)!=scalar(@single_right_file));
}

# open the html template
my $template = HTML::Template->new(filehandle => *DATA);
 
## fill in some parameters
$template->param(max_rd_len => MAX_RD_LEN);

my @loop_lib;
for (my $i=0; $i<=$#left_file; $i++) {
    #print $i, "\t", $left_file[$i], "\t", $right_file[$i];
    
    if (! -e $left_file[$i]) {
        $LOGGER->logdie("Not found paired-end left file ($left_file[$i])");
    }
    if (! -e $right_file[$i]) {
        $LOGGER->logdie("Not found paired-end right file ($right_file[$i])");
    }
    
    my $rf=$left_file[$i];
    $rf=~s/_1\.fastq/_2.fastq/;
    if (basename($right_file[$i]) ne basename($rf)) {
        $LOGGER->logdie("Left (".basename($left_file[$i]).") and Right (".basename($right_file[$i]).") files not match");
    }
    
    my $hr_required = {avg_ins => AVG_INS, rd_len_cutoff => RD_LEN_CUTOFF, pe_left=>abs_path($left_file[$i]), pe_right=>abs_path($right_file[$i] )};
    
    if ($single_left_file[$i]) {
        if (! -e $single_left_file[$i]) {
            $LOGGER->logdie("Not found single-end left file ($single_left_file[$i])");
        }
        my $slf=$left_file[$i];
        $slf=~s/_1\.fastq/_1_singletons.fastq/;
        if (basename($single_left_file[$i]) ne basename($slf)) {
            $LOGGER->logdie("Left (".basename($left_file[$i]).") and Single Left (".basename($slf).") files not match");
        }
        $hr_required->{ 'se_left' } = abs_path($single_left_file[$i]);
        if ($single_right_file[$i]) {
            if (! -e $single_right_file[$i]) {
                $LOGGER->logdie("Not found single-end right file ($single_right_file[$i])");
            }
            my $srf=$right_file[$i];
            $srf=~s/_2\.fastq/_2_singletons.fastq/;
            if (basename($single_right_file[$i]) ne basename($srf)) {
                $LOGGER->logdie("Right (".basename($right_file[$i]).") and Single Right (".basename($srf).") files not match");
            }
            $hr_required->{ 'se_right' } = abs_path($single_right_file[$i]);
        }
    }
    
    push(@loop_lib, $hr_required);
}
$template->param(LIB_LOOP=>\@loop_lib);

# send the obligatory Content-Type and print the template output
print $template->output;


# Subroutines

sub Usage {
    my ($msg) = @_;
	Readonly my $USAGE => <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2019 Universidade Estadual Paulista "Júlio de Mesquita Filho"

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help          Help
        -l      --level         Log level [Default: FATAL]
        -pe1    --pe_left       Input file(s) comma-separated left paired-end reads
        -pe2    --pe_right      Input file(s) comma-separated right paired-end reads
        -se1    --single_left   Input file(s) comma-separated left single-end reads
        -se2    --single_right  Input file(s) comma-separated right single-end reads
        
END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

__DATA__
#https://github.com/alekseyzimin/SOAPdenovo2/blob/master/MANUAL
#maximal read length
max_rd_len=<TMPL_VAR NAME="max_rd_len"><TMPL_LOOP NAME="LIB_LOOP">
[LIB]
#average insert size
avg_ins=<TMPL_VAR NAME="avg_ins">
#if sequence needs to be reversed
reverse_seq=0
#in which part(s) the reads are used
asm_flags=3
#use only first 100 bps of each read
rd_len_cutoff=<TMPL_VAR NAME="rd_len_cutoff">
#in which order the reads are used while scaffolding
rank=1
# cutoff of pair number for a reliable connection (at least 3 for short insert size)
pair_num_cutoff=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
#a pair of fastq file, read 1 file should always be followed by read 2 file
q1=<TMPL_VAR NAME="pe_left">
q2=<TMPL_VAR NAME="pe_right"><TMPL_IF NAME="se_left">
q=<TMPL_VAR NAME="se_left"></TMPL_IF><TMPL_IF NAME="se_right">
q=<TMPL_VAR NAME="se_right"></TMPL_IF></TMPL_LOOP>
