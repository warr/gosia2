Summary: Gosia 2 coulex cross-section code
%define version 20081208.14
%define release 1
Name: gosia2
Version: %{version}
Release: %{release}%{dist}
License: public domain
Vendor: Nigel Warr
Group: Applications/Analysis
Source: /usr/src/redhat/SOURCES/gosia2.tar.gz
Prefix:/usr
BuildRoot: /tmp/package_%{name}-%{version}.%{release}
Requires: gosia_doc

%description
gosia2 is a code for calculating cross-sections for Coulomb excitation
experiments.

%prep
%setup

%build
make

%install
make ROOT="$RPM_BUILD_ROOT" install

%files
%defattr(-,root,root)
/usr/bin/gosia2
/usr/share/man/man1/gosia2.1.gz
