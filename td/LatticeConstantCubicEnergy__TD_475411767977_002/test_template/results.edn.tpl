[{
    "property-id" "tag:staff@noreply.openkim.org,2014-04-15:property/structure-cubic-crystal-npt"
    "instance-id" 1 

    "short-name" {
        "source-value"  ["@<crystal_structure>@"]
    }
    "species" {
        "source-value"  ["@<element>@"]
    }
    "a" {
        "source-value"  @<lattice_constant>@
        "source-unit"   "angstrom"
    }
    "basis-atom-coordinates" {
        "source-value"  @<basis_atoms>@
    }
    "space-group" {
        "source-value"  "@<space_group>@"
    }
    "temperature" {
        "source-value"  0
        "source-unit"   "K"
    }
    "cauchy-stress" {
        "source-value"  [0 0 0 0 0 0]
        "source-unit"   "GPa"
    }
}
{
    "property-id" "tag:staff@noreply.openkim.org,2014-04-15:property/cohesive-potential-energy-cubic-crystal"
    "instance-id" 2  
    "short-name" {
        "source-value"  ["@<crystal_structure>@"]
    }
    "species" {
        "source-value"  ["@<element>@"]
    }
    "a" {
        "source-value"  @<lattice_constant>@
        "source-unit"   "angstrom"
    }
    "basis-atom-coordinates" {
        "source-value"  @<basis_atoms>@
    }
    "space-group" {
        "source-value"  "@<space_group>@"
    }
    "cohesive-potential-energy" {
        "source-value"  @<cohesive_energy>@
        "source-unit"   "eV"
    }
}]
