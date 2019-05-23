#!/bin/perl

# +----------------------------------------------------------------------------+
# |                          Jeffrey M. Hunter                                 |
# |                      jhunter@idevelopment.info                             |
# |                         www.idevelopment.info                              |
# |----------------------------------------------------------------------------|
# |      Copyright (c) 1998-2007 Jeffrey M. Hunter. All rights reserved.       |
# |----------------------------------------------------------------------------|
# | DATABASE : Oracle                                                          |
# | FILE     : dos2unix.pl                                                     |
# | CLASS    : Perl Scripts                                                    |
# | PURPOSE  : Used to convert all end-of-line characters (EOL) of any ASCII   |
# |            file from DOS format to UNIX format. The format for DOS (or     |
# |            Windows) ASCII files use a "crlf" (\r\n) for all EOL            |
# |            characters. UNIX, however, uses only a "lf" (\n) for all EOL    |
# |            characters. When transferring files from DOS to a UNIX system,  |
# |            you will notice ^M characters at the end of each file. This can |
# |            sometimes cause problems with SHELL scripts and configuration   |
# |            scripts. This script can be run on a set of files (or           |
# |            recursively into directories) to replace all "crlf" EOL         |
# |            characters to "lf" only.                                        |
# | NOTE     : As with any code, ensure to test this script in a development   |
# |            environment before attempting to run it in production.          |
# +----------------------------------------------------------------------------+

use File::Find;
use File::Basename;
use Cwd;

# +----------------------------------------------------------------------------+
# | Main Program                                                               |
# +----------------------------------------------------------------------------+

setGlobalVariables();

if ($recursiveSearch) {
    performRecursiveReplace();
} else {
    performReplace();
}


# +----------------------------------------------------------+
# | FUNCTION:    setGlobalVariables                          |
# +----------------------------------------------------------+

sub setGlobalVariables {

    local(@fileList);

    ($scriptName = $0) =~ s#.*/##;
    $recursiveSearch  = 0;
    $backupOld        = 0;
    $find             = "\r";     # find this
    $sub              = undef;    # substitute with this

    foreach $argnum (0 .. $#ARGV) {

        if ($ARGV[$argnum] eq "-r") {
            $recursiveSearch = 1;
        } elsif ($ARGV[$argnum] eq "-b") {
            $backupOld = 1;
        } else {
            push(@fileList, $ARGV[$argnum]);
        }

    }

    @ARGV = @fileList;

    if ($recursiveSearch) {
        print "Performing recursive conversion of all text files.\n";
        print "Ignoring any FILELIST command-line parameters.\n\n";
    } else {
        if (($#fileList + 1) == 0) {
            # print "\nFound 0 arguments for FILELIST.\n";
            usage();
        }
    }

}


# +----------------------------------------------------------+
# | FUNCTION:    usage                                       |
# +----------------------------------------------------------+

sub usage {

    print "\n";
    print "usage: $scriptName [-b] -r || FILELIST\n\n";
    print "       Where:   -b           # backup old file before conversion\n";
    print "                -r           # recurse into directories and convert all text files it finds\n";
    print "                FILELIST     # one or more filenames within the current directory only. '*' = all text files.\n\n";
    exit;

}


# +----------------------------------------------------------+
# | FUNCTION:    performRecursiveReplace                     |
# +----------------------------------------------------------+

sub performRecursiveReplace {

    $From = cwd().'/'.$ARGV[0];
    push(@StartDirectory, $From);
    find(\&displayResults, @StartDirectory);

}


# +----------------------------------------------------------+
# | FUNCTION:    displayResults                              |
# +----------------------------------------------------------+

sub displayResults {

    $fullFileName=$File::Find::name;

    if (-B $fullFileName) {
        ($printfullFileName = $fullFileName) =~ s#$From##;
        if (-d $fullFileName) {
            # printf "Skipping : %-50s   [    directory]\n", $printfullFileName;
        } else {
            printf "Skipping : %-50s   [  binary file]\n", $printfullFileName;
        }
        next;
    }

    ($dev,$ino,$mode,$nlink,$uid,$gid) = stat($fullFileName);

    $backup = $fullFileName . '.bak';
    rename($fullFileName, $backup);

    open(IN, "$backup");
    open(OUT, ">$fullFileName");
    chmod $mode, $fullFileName;

    $counter = 0;

    while($line = <IN>) {
        $counter++;
        $line =~ s/$find/$sub/;
        print OUT $line;
    }

    close(IN);
    close(OUT);

    if ($backupOld) {
        $backupMessage = "Keeping backup";
    } else {
        $backupMessage = "Removing backup";
        unlink $backup;
    }

    ($printfullFileName = $fullFileName) =~ s#$From##;
    # printf "Converted: %-50s   [%7s lines]   [%-15s]\n", $printfullFileName, $counter, $backupMessage;
    printf "Converted: %-50s   [%7s lines]\n", $printfullFileName, $counter;

}


# +----------------------------------------------------------+
# | FUNCTION:    performReplace                              |
# +----------------------------------------------------------+

sub performReplace {

    while (<>) {

        if ($ARGV ne $oldargv) {

            $binaryFile = 0;
            $counter = 0;

            if (-B $ARGV) {
                if (-d $ARGV) {
                    # printf STDOUT "Skipping : %-35s   [    directory]\n", $ARGV;
                } else {
                    printf STDOUT "Skipping : %-35s   [  binary file]\n", $ARGV;
                }
                $oldargv = $ARGV;
                $binaryFile = 1;
                next;
            } else {
                ($dev,$ino,$mode,$nlink,$uid,$gid) = stat($ARGV);
                $backup = $ARGV . '.bak';
                rename($ARGV, $backup);
                open (ARGVOUT, ">$ARGV");
                chmod $mode, $ARGV;
                select(ARGVOUT);
                $oldargv = $ARGV;
            }

        }
        $counter++;
        s/$find/$sub/;

    } continue {

        if ($binaryFile) {

            1;

        } else {

            print;
            if (eof) {
                if ($backupOld) {
                    $backupMessage = "Keeping backup";
                } else {
                    $backupMessage = "Removing backup";
                    unlink $backup;
                }
                # printf STDOUT "Converted: %-35s   [%7s lines]   [%-15s]\n", $oldargv, $counter, $backupMessage;
                printf STDOUT "Converted: %-35s   [%7s lines]\n", $oldargv, $counter;
            }

        }

    }

    select(STDOUT);

}
