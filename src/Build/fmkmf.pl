#!/usr/bin/perl -w
# Copyright (c) The University of Edinburgh
# This is a utility to generate make files 
# for Fortran 90. It was originally in shell script and was re-written
# in perl for greater speed and (hopefully) portability.
# Initial tests suggest speed is 10x better than the sh version.
#
# Format modified by Anna Kelbert, 2005-2009.
#
# A basic makefile entry for bork.f90 would be
# bork.o:bork.f90
# <-tab->$(F90) -c bork.f90
#
# however if bork.f90 contains the line "use gunge" then 
# (A)
# the entry has to be 
# bork.o:bork.f90 garple.o <-- Forces bork to be recompiled if a module it 
# <-tab->$(F90) -c bork.f90                               uses is changed
# where garple.f90 is the program containing the line "module gunge
# (B)  
# The same type of entry has to be done for garple.f90
#
# We also need to generate an entry for the link step. If the main program
# was in baz.f90 then this should be 
# baz:baz.o bork.o.........
# <-tab->$(F90) -o baz baz.o bork.o .....
# The list of object files to be linked should have foo.o in it once 
# and only once for each foo.f90 that was compiled 

use File::Basename;

#-------------------------------------------------
# First check if the luser has any relevent environment vars set
#--------------------------------------------
if( $ENV{FMKMF_F90} ) {
  #print "\# FMKMF_F90 set to $ENV{FMKMF_F90}\n";
  $f90=$ENV{FMKMF_F90};
} 
else {
  #print "\# FMKMF_F90 not set: using f90\n";
  $f90="f90";
}


if( $ENV{FMKMF_SFTAG} ) {
  #print "\# FMKMF_SFTAG set to $ENV{FMKMF_SFTAG}\n";
  $sftag=$ENV{FMKMF_SFTAG};
} 
else {
  #print "\# FMKMF_SFTAG not set: using f90\n";
  $sftag="f90";
}

if( $ENV{FMKMF_SPATH} ) {
  #print "\# FMKMF_SPATH set to $ENV{FMKMF_SPATH}\n";
  $spath=$ENV{FMKMF_SPATH};
} 
else {
  #print "\# FMKMF_SPATH not set: using . \n";
  $spath=".";
}

if( $ENV{FMKMF_OPTIM} ) {
  #print "\# FMKMF_OPTIM set to $ENV{FMKMF_OPTIM}\n";
  $optim=$ENV{FMKMF_OPTIM};
} 
else {
  #print "\# FMKMF_OPTIM not set: using default optimization \n";
  $optim="-O2";
}

if( $ENV{FMKMF_MPIFLAGS} ) {
  #print "\# FMKMF_MPIFLAGS set to $ENV{FMKMF_MPIFLAGS}\n";
  $mpiflags=$ENV{FMKMF_MPIFLAGS};
} 
else {
  #print "\# FMKMF_MPIFLAGS not set: using no link options \n";
  $mpiflags=" ";
}

if( $ENV{FMKMF_LIBRARY_PATH} ) {
  #print "\# FMKMF_LIBRARY_PATH set to $ENV{FMKMF_LIBRARY_PATH}\n";
  $libpath="-L$ENV{FMKMF_LIBRARY_PATH}";
}
else {
  #print "\# FMKMF_LIBRARY_PATH not set: using no library path \n";
  $libpath=" ";
}

if( $ENV{FMKMF_LINKOPTS} ) {
  #print "\# FMKMF_LINKOPTS set to $ENV{FMKMF_LINKOPTS}\n";
  $linkopts=$ENV{FMKMF_LINKOPTS};
} 
else {
  #print "\# FMKMF_LINKOPTS not set: using no link options \n";
  $linkopts=" ";
}

# By default, use the current directory for object files
$linkdir=".";

#------------------------------
# Done with environment variables. Now we need to process commandline args
# These supersede anything supplied via environment variables.
#------------------------------

my ($optiond)=0;
my  $DMPI = 0;

while (@ARGV){

  $arg=shift;
  if($arg =~ /^-win$/){
    $WIN = 1;
    print STDERR "# Creating a Windows style makefile\n";
  }
  if($arg =~ /^-p$/){
    $spath=shift;
    print STDERR "# Using search path from cmd line:\n $spath\n";
  }
  if($arg =~ /^-f90$/){
    $f90=shift;
    print STDERR "# Using compile cmd $f90 from cmd line\n";
  }
  if($arg =~ /^-tag$/){
    $sftag=shift;
    print STDERR "# Using source file tag $sftag from cmd line\n";
  }
  if($arg =~ /^-opt$/){
    $optim=shift;
    print STDERR "# Using compiler optimization options $optim from cmd line\n";
  }
  if($arg =~ /^-mpi$/){
    $mpiflags=shift;
    if($mpiflags =~ /MPI/) {
        $DMPI = 1; # MPI is defined
    }
    print STDERR "# Using compiler MPI flags $mpiflags from cmd line\n";
  }
  if($arg =~ /^-lp$/){
    $libpath=shift;
    print STDERR "# Using Library path $libpath from cmd line\n";
  }
  if($arg =~ /^-l$/){
    $linkopts=shift;
    print STDERR "# Using Link options $linkopts from cmd line\n";
  }
  if($arg =~ /^-d$/){
    $optiond=1;
    print STDERR "# Using debug option (full output on) from cmd line\n";
  }
  if($arg =~ /^-o$/){
      $linkdir=shift;
      print STDERR "# Using $linkdir for object file output directory\n";
  }

}

#-------------------------------------------
# Done processing command line args
#-------------------------------------------


@spath=split(/:/,$spath);


@global_outlines=();
@global_objlist=();
@global_modfiles=();

$mainprogfile=$arg;

# Generate a name for the executable file
$execfile=$mainprogfile;
$execfile=~s/\.${sftag}//;
$execfile=~s|.*/||;

# Output makefile header
print "# Makefile suited for building the $execfile program\n";
print "# Generated using: ./fmkmf.pl [OPTIONS] $mainprogfile > Makefile\n";
print "# with command line options\n";
if($optiond){
  print "# -d (extra debugging output)\n";
}
print "# -p $spath\n";
print "# -f90 $f90 (compiler)\n";
print "# -opt $optim (compiler optimisation)\n";
print "# -lp $libpath (linking options: path to libraries)\n";
print "# -l $linkopts (linking options)\n";
print "# -o $linkdir (output directory for object files)\n\n";

if($optiond){
  print STDERR "# Main program is $mainprogfile \n" ;
}
# this subroutine (def below) does most of the work.
process_fsource($mainprogfile); 

# set some makefile . 

print "\n# ------------------Macro-Defs---------------------\n";

print "include Build/Makefile.local\n";
print "OBJDIR = ../../objs/MT3D_OO\n";
print "F90 = $f90\n";
print "FFLAGS = $optim\n";
print "MPIFLAGS = $mpiflags\n";
if($WIN) {
    print "MODULE = \n";
} elsif($f90 =~ /^g95$/){
    print "MODULE = -fmod=\$(OBJDIR)\n";
} elsif($f90 =~ /^gfortran$/){
    print "MODULE = -J \$(OBJDIR)\n";
} else {
    print "MODULE = -J \$(OBJDIR)\n";
}
if($libpath !~ /^(\s*)$/){
    print "LIBS_PATH = -L$libpath\n";
} else {
    print "LIBS_PATH = \n";    
}
print "LIBS = $linkopts\n";


print "\n# -------------------End-macro-Defs---------------------------\n";

print "OBJ = @global_objlist \n\n";


print "\nall: $execfile \n";

# Generate makefile entry for the Link step
print "\n# Here is the link step \n";

if($WIN) {
    print "$execfile: \$(OBJDIR) \$(OBJ) \n";
    print "\t \$(F90) /link \$(OUTDIR)/$execfile \$(OBJ) \$(LIBS_PATH) \$(LIBS)\n";
} else {
    print "$execfile: \$(OBJDIR) \$(OBJ) \n";
    print "\t \$(F90) -o \$(OUTDIR)/$execfile \$(OBJ) \$(LIBS_PATH) \$(LIBS)\n";
}
# print "\trm -f *.mod \n";

print "\n# Here are the compile steps \n\n";

print "\$(OBJDIR): \n";
print "\tmkdir -p \$(OBJDIR)\n";

print STDOUT @global_outlines;

# Add an entry for make clean at the end of the make file.  this
# removes most of the garbage left around by most of the Fortran 90
# compilers I have tried.

print "\n# Type \" make clean \" to get rid of all object and module files \n";

print "clean: \n";
print "\tcd \$(OBJDIR); \\\n";
print "\trm -f *~ *.o *.obj *.mod *.d *.s00 *.dbg *.stackdump \\\n";
print "\t`find . -mindepth 1 -name \"*~\"` \n\n";
print "cleanall: clean \n";
print "\trm -f \$(OUTDIR)/$execfile \n\n";
print "src: clean \n";
print "\ttar cvfz \$(ARCHIVE).tgz * \n";
print "\n  \n";

# End of main program 

##############################################
# Here is the subroutine that generates the compile entries in the makefile
# These end up in the global array @global_outlines. The magic part is 
# that this subroutine calls itself recursively.
##############################################
sub process_fsource {
  
  my $mainprogfile=$_[0];
  if($optiond){
    print STDERR "# process_fsource called with arg $mainprogfile \n";
  }
  open( MAINPROG, $mainprogfile) or 
    die "Can't find main program file $mainprogfile: $!> \n";
  
  # Read through Fortran source looking for USE statements
  # There should be nothing but whitespace before the USE. Sloppily, 
  # we allow tabs, although the standard (IIRC) does not
  # An important exception, helpful for our purposes:
  # if the used module has MPI in the name, omit it unless MPI is defined.
  # In the code, this is accomplished with #ifdef MPI statements.
  my @modulelist=();
  my @includelist=();
  while ($line=<MAINPROG>) { 
    if($line =~ /^[ \t]*use (\w+)/i ) { # line matches regexp between / /
      my $modulefile = $1;
      if($optiond){
        print STDERR "# $mainprogfile Uses Module $modulefile\n";
      }
      if($modulefile =~ /MPI/) { # if MPI is found in module name, skip unless MPI is defined
          next unless ($DMPI);
      }
      @modulelist=(@modulelist,$modulefile);
    } elsif  ($line =~ /^[ \t]*\#?include "(\S+)"/i ){
      my $includefile = $1;
      if($optiond){
        print STDERR "# $mainprogfile Includes $includefile\n";
      }
      if($includefile =~ /MPI/) { # if MPI is found in module name, skip unless MPI is defined
          next unless ($DMPI);
      }
      my $mainprogpath = dirname($mainprogfile); 
      @includelist=(@includelist,"$mainprogpath/$includefile");        
    }
  }
  
  close(MAINPROG);

  if($optiond){
    print STDERR "# Full list of modules in $mainprogfile: @modulelist \n";
    if(@includelist > 0) {
        print STDERR "# Full list of includes in $mainprogfile: @includelist \n";
    }
  }
  # Find which file each module is in.
  
  
  
 my @modfiles=();
 MODLOOP:foreach $module (@modulelist){
    foreach $directory (@spath){
      # print "# Looking in directory $directory\n";
      opendir( DIRHANDLE, $directory) or die 
    "Can't open directory $directory : $!> \n";
      @sourcefiles=grep /\.${sftag}\Z/, sort(readdir(DIRHANDLE));
    foreach $sourcefile (@sourcefiles){
      $pathsourcefile="$directory/$sourcefile";
      #print "\# Checking $pathsourcefile\n";
      open( SOURCEFILE, "$pathsourcefile") or 
    die "Can't find source file $pathsourcefile: $!> \n";
      while ($line=<SOURCEFILE>){
    if($line =~ /^ *module (\w+)/i ){
      if($1 =~ /^$module$/i){
        if($optiond){
          print STDERR "# Uses $module which is in $pathsourcefile\n";
        }
        @modfiles=(@modfiles,$pathsourcefile);
        
        if(grep (/$pathsourcefile/,@global_modfiles )){
          if($optiond){
        print STDERR "# $pathsourcefile already in list\n";
          }
        }
        else {
          @global_modfiles=(@global_modfiles,$pathsourcefile);
          process_fsource($pathsourcefile);

        }
        # We found this module -- go on to the next one
        close (SOURCEFILE);
        next MODLOOP;        
      }
    }
      }
      close( SOURCEFILE );
    }
  }
  # exhausted source files
  print STDERR "Couldn't find source file for module $module\n";
}

if($WIN) {
    $mainprogfile=~s/\//\\/g;
    # name of file we want to make
    $objfile=$mainprogfile;
    # replace source file name with .o
    $objfile=~s/\.${sftag}/\.obj/;
    # strip path so object files go in current dir
    $objfile=~s|.*\\||;
} else {
    # name of file we want to make
    $objfile=$mainprogfile;
    # replace source file name with .o
    $objfile=~s/\.${sftag}/\.o/;
    # strip path so object files go in current dir
    $objfile=~s|.*/||;
}

# now add the user-defined path to the object files
$objfile="\$(OBJDIR)/$objfile";
@global_objlist=(@global_objlist,$objfile);
# list of dependencies
@objlist=();
foreach  $mf (@modfiles) { 
  $obj=$mf;
  # replace source file name with .o
  if($WIN) {
      $obj=~s/\.${sftag}/\.obj/;
  } else {
      $obj=~s/\.${sftag}/\.o/;
  }
  # strip path so object files go in current dir
  $obj=~s|.*/||;
  # now add the user-defined path to the object file
  $obj="\$(OBJDIR)/$obj";
  @objlist=(@objlist,$obj);
}

@global_outlines=(@global_outlines,"\n$objfile:$mainprogfile @objlist @includelist\n");
if($WIN){
    @global_outlines=(@global_outlines,"\t \$(F90) -c \$(MODULE) \$(FFLAGS) \$(MPIFLAGS) \"$mainprogfile\" /link $objfile\n");
} else {
    @global_outlines=(@global_outlines,"\t \$(F90) -c \$(MODULE) \$(FFLAGS) \$(MPIFLAGS) $mainprogfile -o $objfile\n");
}
#if(@includelist > 0) {
#    @global_outlines=(@global_outlines,"\n$mainprogfile: @includelist \n");
#}

}
