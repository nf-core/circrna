# barriers.spec
#
# Copyright (c) 2001 Ivo Hofacker, Peter Stadler, ivo@tbi.univie.ac.at
#
%define name ViennaRNA
%define version 1.6
%define release 1
%define manifest %{_builddir}/%{name}-%{version}-%{release}.manifest

# required items
Name: %{name}
Version: %{version}
Release: %{release}
Copyright: GPL
Group: Application/Misc

# optional items
Vendor: Ivo Hofacker, Peter Stadler, TBI
#Distribution:
#Icon:
URL: http://www.tbi.univie.ac.at/~ivo/RNA/
Packager: Ivo Hofacker, ivo@tbi.univie.ac.at

# source + patches
Source: %{name}-%{version}.tar.gz
#Source1:
#Patch:
#Patch1:

# RPM info
Provides: libRNA.a
#Requires:
#Conflicts:
#Prereq:

#Prefix: /usr
BuildRoot: /var/tmp/%{name}-%{version}

Summary: RNA secondary structure prediction and analysis.

%description
The ViennaRNA package consists of a library and several standalone
programs for RNA secondary structure analysis. It includes algorithms
for predicting optimal and suboptimal secondary structures, base pair
probabilities and partition functions, for comparing secondary
structures, and the design of RNA sequences with a desired structure.
 

%prep
%setup -q
#%patch0 -p1

%build
%configure
make

%install
rm -rf $RPM_BUILD_ROOT
mkdir -p $RPM_BUILD_ROOT
%makeinstall

# __os_install_post is implicitly expanded after the
# %install section... do it now, and then disable it,
# so all work is done before building manifest.

%{?__os_install_post}
%define __os_install_post %{nil}

# build the file list automagically into %{manifest}

cd $RPM_BUILD_ROOT
rm -f %{manifest}
find . -type d \
        | sed '1,2d;s,^\.,\%attr(-\,root\,root) \%dir ,' >> %{manifest}
find . -type f \
        | sed 's,^\.,\%attr(-\,root\,root) ,' >> %{manifest}
find . -type l \
        | sed 's,^\.,\%attr(-\,root\,root) ,' >> %{manifest}

#%pre
#%post
#%preun
#%postun

%clean
rm -f %{manifest}
rm -rf $RPM_BUILD_ROOT

%files -f %{manifest}
%defattr(-,root,root)
#%doc README
#%docdir
#%config

%changelog
