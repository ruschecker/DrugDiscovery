#!/usr/bin/perl -w

my $ARG="";
$receptor_gro=$ARGV[0];
print "MDRUN PREP STARTED\n";

report_date();
running_check();

#Set Variables
#Project Directory
$nsteps = 2500000;
$steps=50000;
$start=0;
$end = $start + $steps;
$nstxout=($steps * 0.10);
$nstvout=($steps * 0.10);
$nstxtcout=($steps * 0.10);
$rsc_fpops_est = ($steps * 319500000);
$rsc_fpops_bound = ($rsc_fpops_est * 1000);
$rsc_memory_bound = 500000000;
$delay_bound = ($steps * 7);
$PROJECT="/home/boincadm/projects/DrugDiscovery";
$OPS=$PROJECT."/html/ops";
$results=$PROJECT."/sample_results";


$time=`date '+%s%N'`;
chomp ($time);

$job_file = "
<job_desc>
<task>
<application>mdrun.exe</application>
<command_line> -v -x -c -o -e -cpi md.cpt -cpt 1 -deffnm md</command_line>
<weight>100</weight>
<fraction_done_filename>progress.txt</fraction_done_filename>
</task>
<task>
<application>7za.exe</application>
<weight>1</weight>
<command_line> a -y out.7z Ligand.acpypi.7z *.log *.tpr *.gro *.mdp *.edr *.cpt job.xml *.xtc *.top *.itp</command_line>
</task>
</job_desc>
";

$make_ndx = "0 q ";
$protein_lig = "protein LIG ";
$g_energy = "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ";

#use the Cwd methods
use Cwd;
#save the original directory path
my $orig_dir = Cwd::abs_path;
#clean the path string
chomp($OPS);
#change to new directory
chdir($OPS);

#change to original 
chdir($orig_dir);

#remove old assignment files

#Sample Results Directory
$results=$PROJECT."/sample_results";
$tmp=$results."/tmp";

if (-d $tmp) {
  $ARG="rm -rf ".$tmp;
  system($ARG);
}

chdir($results);
$ARG="mkdir tmp";
system($ARG);

#Change to cp for asgn
$ARG="cp ".$results."/autodock*.pdbqt tmp/";
system($ARG);

chdir($tmp);
$pwd=`pwd`;
print $pwd."\n";

$ARG="/usr/local/bin/pythonsh ".$PROJECT."/bin/vina_screen_get_top.py 100 >".$tmp."/top_autodock.txt";
print $ARG."\n";
system($ARG);

$ARG="sed -i '/setting\ PYTHONHOME\ environment/d' ".$tmp."/top_autodock.txt";
print $ARG."\n";
system($ARG);

$top_scores = $tmp."/top_autodock.txt";
open (TOP_SCORES,$top_scores) or die ("Error trying to open top scores.\n");
# process top scores

do {$new_score = <TOP_SCORES>;

  print $new_score."\n";

  report_date();

  #generate seed
  #  use strict;
  #  use warnings;
  my $range = 1000000;
  my $seed = int(rand($range));

  $md_param="integrator \= md
  nsteps \= ".$steps."
  dt \= 0.002
  nstvout \= ".$nstvout."
  nstlog \= 500
  nstenergy \= 250
  nstxtcout \= ".$nstxtcout."
  nstxout \= ".$nstxout."
  xtc_grps \= Protein SOL
  energygrps \= Protein  SOL
  constraints \= all-bonds
  nstcomm \= 1
  ns_type \= grid
  rlist \= 1.2
  rcoulomb \= 1.1
  rvdw \= 1.0
  vdwtype \= shift
  rvdw-switch \= 0.9
  coulombtype \= PME-Switch
  Tcoupl \= v-rescale
  tau_t \= 0.1 0.1
  tc-grps \= protein non-protein
  ref_t \= 300 300
  Pcoupl \= parrinello-rahman
  PcOupltype \= isotropic
  tau_p \= 0.5
  compressibility \= 4.5e-5
  ref_p \= 1.0
  gen_vel \= yes
  lincs-iter \= 2
  DispCorr \= EnerPres
  optimize_fft \= yes
  gen_seed = ".$seed;


  $df=`df | head -n 3 | tail -n 1 | awk \'{ print substr( \$0, length(\$0) - 4, 2 ) }\'`;
  chomp($df);
  $df=$df-1;
  if ($df >70) {
          print "over 70% full \n";
          sleep(5200);
  }

  chomp($new_score);
  #if not at the end of file
  if ($new_score ne "EOF")
  {
    $time=`date '+%s%N'`;
    chomp ($time);

    #create a new score directory
    $ARG="mkdir ".$new_score.".tmp";
    print $ARG."\n";
    system($ARG);
    $ARG="cp ".$new_score." ".$new_score.".tmp/ligand.pdbqt";
    print $ARG."\n";
    system($ARG);
    #use the Cwd methods
    use Cwd;
    #save the original directory path
    my $orig_dir = Cwd::abs_path;
    #clean the path string
    chomp($new_score);
    #change to new directory
    chdir($new_score.".tmp");

    $bt=substr($new_score, rindex($new_score,"_bt_")+4, rindex($new_score,"_lig_") - rindex($new_score,"_bt_")-4);
    $lig=substr($new_score, rindex($new_score,"_lig_")+5, rindex($new_score,"_ts_") - rindex($new_score,"_lig_")-5);
    $ts=substr($new_score, rindex($new_score,"_ts_")+4);
    $wu_name="nsteps_".$nsteps."_int_".$start."_-_".$end."_bt_".$bt."_lig_".$lig."_ts_".$ts;

    $ENV{'PATH'} = '/usr/local/gromacs/bin/:/usr/kerberos/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/antechamber-1.27:/usr/local/antechamber-1.27/exe:/home/boincadm/bin';
    $ENV{'ACHOME'} = '/usr/local/antechamber-1.27';
    $ENV{'AMBERHOME'} = '/usr/local/antechamber-1.27/exe';


    $file=$new_score.".tmp/ligand.pdbqt";
    #Calculate file size

    $filesize = -s $file;
    #Error handle: If file size is not greater than 10 bytes, skip this file

    #Convert PDBQT format to PDB
    $ARG="/usr/local/bin/pythonsh /usr/local/MGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/pdbqt_to_pdb.py -f ligand.pdbqt -o ligand.pdb";        
    print $ARG."\n";
    system($ARG);
    #Babel adds protons
    $ARG="babel -ipdb ligand.pdb -opdb ligand.pdb -h";
    print $ARG."\n";
    system($ARG);
    $ARG="sed -i \'s\/<0>\/LIG\/g\' ligand.pdb";
    print $ARG."\n";
    system($ARG);


          #Output only LIG residues
  $ARG="grep \"LIG\" ligand.pdb > Ligand.pdb";
  print $ARG."\n";
  system($ARG);

          #Acpypi preps the Ligand file for GROMACS
  $ARG="/usr/local/bin/acpypi -f -d -n 0 -i Ligand.pdb -s 120";
  print $ARG."\n";
  system($ARG);

  $filename=$new_score."/Ligand.acpypi/Ligand_GMX.itp";
  if (-e $filename) {

  #Copy the Protein file previously preped with Amber
  $ARG="cp ".$PROJECT."/biotarget/".$receptor_gro." .";
  print $ARG."\n";
  system($ARG);


  #Prep the Protein2.pb useing pdb2gmx
  $ARG="pdb2gmx -ff amber99sb -f ".$receptor_gro." -o Protein2.pdb -p Protein.top -water spce -ignh -missing";
  print $ARG."\n";
  system($ARG);

exit;

  #Extract the ATOM Residues from Protein2.pdb and Ligand_NEW.pdb and redirect to Complex.pdb
  $ARG="grep -h ATOM Protein2.pdb Ligand.acpypi/Ligand_NEW.pdb >| Complex.pdb";
  print $ARG."\n";
  system($ARG);

  #Copy the Ligand_GMX.itp to Ligand.itp
  $ARG="cp Ligand.acpypi/Ligand_GMX.itp Ligand.itp";
  print $ARG."\n";
  system($ARG);

  #Copy the Protein topology file to Complex.top
  $ARG="cp Protein.top Complex.top";
  print $ARG."\n";
  system($ARG);

  #Copy the Complex.top and insert #include "Ligand.itp" into new file Complex2.top
  $ARG="cat Complex.top | sed '/#include\ \\\"ffamber99sb.itp\\\"/a \#include \"Ligand.itp\"' >| Complex2.top";
  print $ARG."\n";
  system($ARG);

  #Append Ligand   1 to Complex2.top
  $ARG="echo \"Ligand   1\" >> Complex2.top";
  print $ARG."\n";
  system($ARG);

  #Rename the Complex2.top to Complex.top
  $ARG="mv Complex2.top Complex.top";
  print $ARG."\n";
  system($ARG);

  #Copy the em.mdp parameter file to current working directory
  $ARG="cp ".$PROJECT."/bin/em.mdp .";
  print $ARG."\n";
  system($ARG);

  $ARG="cp ".$PROJECT."/bin/em_r.mdp .";
  print $ARG."\n";
  system($ARG);

  $ARG="cp ".$PROJECT."/bin/pr.mdp .";
  print $ARG."\n";
  system($ARG);

  #Run editconf on the Complex.pdb
  $ARG="editconf -bt triclinic -f Complex.pdb -o Complex.pdb -d 1.0";
  print $ARG."\n";
  system($ARG);

  #Run genbox on Complex.pdb
  $ARG="genbox -cp Complex.pdb -cs ffamber_tip3p.gro -o Complex_b4ion.pdb -p Complex.top";
  print $ARG."\n";
  system($ARG);

  #Run grompp
  $ARG="grompp -f em_r.mdp -c Complex_b4ion.pdb -p Complex.top -o Complex_b4ion.tpr";
  print $ARG."\n";
  system($ARG);

  #copy the Complex.top to Complex_ion.top
  $ARG="cp Complex.top Complex_ion.top";
  print $ARG."\n";
  system($ARG);

  #Run genion
  $ARG="echo 13| genion -s Complex_b4ion.tpr -o Complex_b4em.pdb -neutral -conc 0.15 -p Complex_ion.top -norandom";
  print $ARG."\n";
  system($ARG);

  sleep(1);

  #rename the Complex_ion.top Complex.top
  $ARG="mv Complex_ion.top Complex.top";
  print $ARG."\n";
  system($ARG);
  #Run the grompp
  $ARG="grompp -f em_r.mdp -c Complex_b4em.pdb -p Complex.top -o em_r.tpr";
  print $ARG."\n";
  system($ARG);

  my $timeout = 180;
  my $pid = fork;

  if ( defined $pid ) {
  if ( $pid ) {
  # this is the parent process
  local $SIG{ALRM} = sub { die "TIMEOUT" };
  alarm 3600;
  # wait until child returns or timeout occurs
  eval {
          waitpid( $pid, 0 );
  };
  alarm 0;

  if ( $@ && $@ =~ m/TIMEOUT/ ) {
          # timeout, kill the child process
          kill 9, $pid;
  }
  }
  else {
  # this is the child process
  # this call will never return. Note the use of exec instead of run_system
  exec "mdrun -v -deffnm em_r -maxh 1";
  }
  }
  else {
          die "Could not fork.";
  }
  $ARG="grompp -f em.mdp -c em_r.gro -p Complex.top -o em.tpr";
  system($ARG);

  my $pid = fork;

  if ( defined $pid ) {
  if ( $pid ) {
  # this is the parent process
  local $SIG{ALRM} = sub { die "TIMEOUT" };
  alarm 36000;
  # wait until child returns or timeout occurs
  eval {
          waitpid( $pid, 0 );
  };
  alarm 0;

  if ( $@ && $@ =~ m/TIMEOUT/ ) {
          # timeout, kill the child process
          kill 9, $pid;
  }
  }
  else {
  # this is the child process
  # this call will never return. Note the use of exec instead of run_system
  exec "mdrun -v -deffnm em -maxh 1";
  }
  }
  else {
          die "Could not fork.";
  }

  $ARG="grompp -f pr.mdp -c em.gro -p Complex.top -o pr.tpr";
  system($ARG);

  my $pid = fork;

  if ( defined $pid ) {
  if ( $pid ) {
  # this is the parent process
  local $SIG{ALRM} = sub { die "TIMEOUT" };
  alarm 3600;
  # wait until child returns or timeout occurs
  eval {
          waitpid( $pid, 0 );
  };
  alarm 0;

  if ( $@ && $@ =~ m/TIMEOUT/ ) {
          # timeout, kill the child process
          kill 9, $pid;
  }
  }
  else {
  # this is the child process
  # this call will never return. Note the use of exec instead of run_system
  exec "mdrun -v -s pr.tpr -o pr.trr -c b4md.gro -g pr.log -e pr.edr";
  }
  }
  else {
          die "Could not fork.";
  }


  open (MYFILE, ">>".$new_score."/md.mdp");
  close (MYFILE);

  $ARG="grompp -f md.mdp -c b4md.gro -p Complex.top -o md.tpr";
  system($ARG);

  open (MYFILE, ">>".$new_score."/md_job.xml");
  close (MYFILE);

  $ARG="7za a -r Ligand.acpypi.7z Ligand.acpypi";
  system($ARG);


  chdir($PROJECT);

  #copy the md.tpr used by mdrun workunit. Copy to the download directory using the dir_hier_path which sets the proper fanout directory
  $ARG="cp ".$new_score."/md.tpr \`bin/dir_hier_path md_".$wu_name.".tpr\`";
  system($ARG);

  $ARG="cp ".$new_score."/md.mdp \`bin/dir_hier_path md_".$wu_name.".mdp\`";
  system($ARG);

  $ARG="cp ".$new_score."/md_job.xml \`bin/dir_hier_path md_job_".$wu_name.".txt\`";
  system($ARG);

  $ARG="cp ".$new_score."/Complex.top \`bin/dir_hier_path complex".$wu_name.".txt\`";
  system($ARG);

  $ARG="cp ".$new_score."/Ligand.acpypi.7z \`bin/dir_hier_path Ligand.acpypi_".$wu_name.".7z\`";
  system($ARG);

  $ARG="cp ".$new_score."/Ligand.itp \`bin/dir_hier_path Ligand_".$wu_name.".itp\`";
  system($ARG);

  $ARG="chmod -R a+r ".$PROJECT."/download";
  system($ARG1);


  $ARG=$PROJECT."/bin/create_work -appname mdrun -wu_name ".$wu_name." -wu_template ".$PROJECT."/templates/mdrun_wu -result_template templates/mdrun_result --rsc_fpops_est ".$rsc_fpops_est." --rsc_fpops_bound ".$rsc_fpops_bound." --delay_bound ".$delay_bound." md_".$wu_name.".tpr md_".$wu_name.".mdp md_job_".$wu_name.".txt complex".$wu_name.".txt Ligand.acpypi_".$wu_name.".7z Ligand_".$wu_name.".itp";

  system($ARG);
  system("touch ".$PROJECT."/reread_db");

  report_date();

  }

  chdir($orig_dir);
  $ARG="rm -rf ".$new_score." ".$new_score.".7z";
  system($ARG);
  }
  } while ($new_score ne "EOF");
  close (TOP_SCORE);

  report_date();

  chdir($results);

  $ARG="mv ".$tmp."/*.7z ".$PROJECT."/backup";
  system($ARG);

  system("rm -rf ".$tmp);

  sub running_check
  {
  $check = system("ps ax | grep mdrun.pl | wc -l");
  chomp($check);
  print "Checking if script is running\n";
  if ($check > 10) {
          print "Already running\n";
          exit 0;
  }
  }

  sub report_date {
  $ARG="date";
  system($ARG);

  }

  sub system 
  {
  system($ARG);
  if ( $? == -1 )
  {
  print "command failed: $!\n";
  }
  else
  {
  printf $ARG.": command exited with value %d", $? >> 8;
  }
  }
