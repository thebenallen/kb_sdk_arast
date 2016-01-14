#! /usr/bin/env perl

use strict;

use Carp;
use Config::Simple;
use Cwd 'abs_path';
use Data::Dumper;
use File::Basename;
use JSON 'decode_json';

my $repo_dir = abs_path(dirname(__FILE__)."/..");

my $usage = "Usage: $0 plugin assembly.dir\n\n";
my $plugin = shift @ARGV or die $usage;

my $asm_dir = shift @ARGV || "$repo_dir/../../assembly";
my $info = plugin_info($plugin, $asm_dir);
my $desc = $info->{desc};
my $params = $info->{params};
my $param = '('.join(' ', map { "$_=$params->{$_}" } keys %$params).')';
# print STDERR '$param = '. Dumper($param);
# print STDERR '$info = '. Dumper($info);

my $spec_dir = "$repo_dir/ui/narrative/methods/run_$plugin";
verify_dir($spec_dir);

print STDERR "Create $spec_dir/spec.json\n";
run("tpage --define plugin=$plugin $repo_dir/tools/spec.json.tt > $spec_dir/spec.json");

print STDERR "Create $spec_dir/display.yaml\n";
run("tpage --define plugin=$plugin --define description='$desc' --define params='$param' $repo_dir/tools/display.yaml.tt > $spec_dir/display.yaml");
append_pubs($spec_dir, $info->{refs});

print STDERR "Append $repo_dir/lib/AssemblyRAST/AssemblyRASTImpl.py\n";
run("tpage --define plugin=$plugin $repo_dir/tools/AssemblyRASTImpl.py.tt >> $repo_dir/lib/AssemblyRAST/AssemblyRASTImpl.py") if !`grep run_$plugin $repo_dir/lib/AssemblyRAST/AssemblyRASTImpl.py`;

print STDERR "Update $repo_dir/AssemblyRAST.spec\n";
update_module_spec($plugin, "$repo_dir/AssemblyRAST.spec");

sub update_module_spec {
    my ($plugin, $file) = @_;
    return if `grep run_$plugin $file`;
    my $spec = `cat $file`;
    $spec =~ s/\s*};\s*$//;
    $spec .= "\n\n";
    $spec .= <<"End_of_funcdef";
    funcdef run_$plugin(AssemblyParams params) returns (AssemblyOutput output)
        authentication required;

};
End_of_funcdef
    open(SPEC, ">$file") or die "Could not open $file";
    print SPEC $spec;
    close(SPEC);
}

sub append_pubs {
    my ($spec_dir, $refs) = @_;
    my $yaml = "$spec_dir/display.yaml";
    my $pubs = refs_to_yaml($info->{refs});
    open(YAML, ">>$yaml") or die "Could not open $yaml";
    print YAML "\npublications :\n";
    print YAML $pubs;
    close(YAML);
}

sub refs_to_yaml {
    my ($refs) = @_;
    my @yaml;
    for my $ref (@$refs) {
        push @yaml, "    -";
        if ($ref =~ /^doi:(\S+)/) {
            my $citation = doi_to_citation($1);
            push @yaml, "        display-text : |";
            push @yaml, "            $citation";
        } elsif ($ref =~ /^http/) {
            push @yaml, "        link: $ref";
        } else {
            push @yaml, "        display-text : |";
            push @yaml, "            $ref";
        }
    }
    return join("\n", @yaml) . "\n";
}

sub doi_to_hash {
    my ($doi) = @_;
    $doi =~ s/^doi://;
    my $json = `curl -s -L -H "Accept: application/json" "http://dx.doi.org/$doi"`;
    return decode_json($json);
}

sub doi_to_citation {
    my ($doi) = @_;
    my $hash = doi_to_hash($doi);
    # print STDERR '$hash = '. Dumper($hash);

    # Example: Zerbino, D. R., & Birney, E. (2008). Velvet: algorithms for de novo short read assembly using de Bruijn graphs. Genome research, 18(5), 821-829.
    my $author = join(", ", map { join(', ', $_->{family}, $_->{given}) } @{$hash->{author}});
    my $year = $hash->{issued}->{'date-parts'}->[0]->[0];
    my $title = $hash->{title};
    my $journal = $hash->{'container-title'};
    my $volume = $hash->{volume};
    my $issue = $hash->{issue};
    my $page = $hash->{page};
    $journal .= ", $volume($issue), $page" if $volume && $issue && $page;
    my $citation = join(" ", $author, "($year)", $title.'.', $journal) if $title;
    $citation .= ", doi: $doi";
    return $citation;
}

sub plugin_info {
    my ($plugin, $asm_dir) = @_;
    -d $asm_dir or die "Invalid assembly repo path: '$asm_dir'\n\n$usage";
    my $file = "$asm_dir/lib/assembly/plugins/$plugin.asm-plugin";
    -s $file or die "Plugin implementation not found: $file\n\n$usage";
    my $cfg = new Config::Simple($file);
    my %info;
    $info{desc} = $cfg->param("Documentation.Description");
    $info{params} = $cfg->get_block('Parameters'); # hash
    $info{refs} = [ split(/[, ]+/, $cfg->param("Documentation.References")) ];
    wantarray ? %info : \%info;
}

sub verify_dir {
    my ($dir) = @_;
    run("mkdir -p $dir");
}

sub run { system(@_) == 0 or confess("FAILED: ". join(" ", @_)); }
