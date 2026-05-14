#!/usr/bin/env julia

using TOML
using Random
using Printf
using ArgParse

# Types
mutable struct Particle
  mol::String
  atm::String
  mid::Int64
  aid::Int64
  x::Float64
  y::Float64
  z::Float64
end

# Force field and properties
ISOLF = Dict(
  "dlpa" => Dict("beads" => [       "PHA", "MID", "DL1", "DL2"       ], "bonds" => [        0.4685, 0.4745, 0.5645        ], "phosp" => 1, "apl" => 0.5468),
  "dlpc" => Dict("beads" => ["CHO", "PHO", "MID", "DL1", "DL2"       ], "bonds" => [0.4025, 0.4715, 0.4675, 0.5505        ], "phosp" => 2, "apl" => 0.6199),
  "dlpe" => Dict("beads" => ["ETH", "PHO", "MID", "DL1", "DL2"       ], "bonds" => [0.3565, 0.4735, 0.4755, 0.5665        ], "phosp" => 2, "apl" => 0.5316),
  "dlpg" => Dict("beads" => ["GLY", "PHO", "MID", "DL1", "DL2"       ], "bonds" => [0.3815, 0.4695, 0.4665, 0.5455        ], "phosp" => 2, "apl" => 0.6489),
  "dlps" => Dict("beads" => ["SER", "PHO", "MID", "DL1", "DL2"       ], "bonds" => [0.4425, 0.4695, 0.4735, 0.5635        ], "phosp" => 2, "apl" => 0.5448),

  "dmpa" => Dict("beads" => [       "PHA", "MID", "DM1", "DM2"       ], "bonds" => [        0.4705, 0.5515, 0.7245        ], "phosp" => 1, "apl" => 0.4379),
  "dmpc" => Dict("beads" => ["CHO", "PHO", "MID", "DM1", "DM2"       ], "bonds" => [0.4035, 0.4715, 0.5235, 0.6635        ], "phosp" => 2, "apl" => 0.5965),
  "dmpe" => Dict("beads" => ["ETH", "PHO", "MID", "DM1", "DM2"       ], "bonds" => [0.3555, 0.4785, 0.5545, 0.7565        ], "phosp" => 2, "apl" => 0.4314),
  "dmpg" => Dict("beads" => ["GLY", "PHO", "MID", "DM1", "DM2"       ], "bonds" => [0.3815, 0.4715, 0.5195, 0.6515        ], "phosp" => 2, "apl" => 0.6238),
  "dmps" => Dict("beads" => ["SER", "PHO", "MID", "DM1", "DM2"       ], "bonds" => [0.4435, 0.4725, 0.5415, 0.6995        ], "phosp" => 2, "apl" => 0.4597),

  "dppa" => Dict("beads" => [       "PHA", "MID", "DP1", "DP2", "DP3"], "bonds" => [        0.4725, 0.4995, 0.6345, 0.5745], "phosp" => 1, "apl" => 0.4407),
  "dppc" => Dict("beads" => ["CHO", "PHO", "MID", "DP1", "DP2", "DP3"], "bonds" => [0.4025, 0.4755, 0.4895, 0.6315, 0.5725], "phosp" => 2, "apl" => 0.5185),
  "dppe" => Dict("beads" => ["ETH", "PHO", "MID", "DP1", "DP2", "DP3"], "bonds" => [0.3555, 0.4795, 0.4985, 0.6345, 0.5745], "phosp" => 2, "apl" => 0.4363),
  "dppg" => Dict("beads" => ["GLY", "PHO", "MID", "DP1", "DP2", "DP3"], "bonds" => [0.3805, 0.4735, 0.4915, 0.6315, 0.5735], "phosp" => 2, "apl" => 0.5018),
  "dpps" => Dict("beads" => ["SER", "PHO", "MID", "DP1", "DP2", "DP3"], "bonds" => [0.4425, 0.4745, 0.4975, 0.6335, 0.5745], "phosp" => 2, "apl" => 0.4613),

  "dopa" => Dict("beads" => [       "PHA", "MID", "DO1", "DO2", "DO3"], "bonds" => [        0.4665, 0.4645, 0.5085, 0.5575], "phosp" => 1, "apl" => 0.6359),
  "dopc" => Dict("beads" => ["CHO", "PHO", "MID", "DO1", "DO2", "DO3"], "bonds" => [0.4025, 0.4705, 0.4615, 0.5035, 0.5475], "phosp" => 2, "apl" => 0.6793),
  "dope" => Dict("beads" => ["ETH", "PHO", "MID", "DO1", "DO2", "DO3"], "bonds" => [0.3565, 0.4695, 0.4645, 0.5155, 0.5605], "phosp" => 2, "apl" => 0.6179),
  "dopg" => Dict("beads" => ["GLY", "PHO", "MID", "DO1", "DO2", "DO3"], "bonds" => [0.3805, 0.4695, 0.4575, 0.4975, 0.5375], "phosp" => 2, "apl" => 0.7090),
  "dops" => Dict("beads" => ["SER", "PHO", "MID", "DO1", "DO2", "DO3"], "bonds" => [0.4425, 0.4665, 0.4655, 0.5095, 0.5575], "phosp" => 2, "apl" => 0.6330),

  "dspa" => Dict("beads" => [       "PHA", "MID", "DS1", "DS2", "DS3"], "bonds" => [        0.4725, 0.4995, 0.6345, 0.7005], "phosp" => 1, "apl" => 0.4492),
  "dspc" => Dict("beads" => ["CHO", "PHO", "MID", "DS1", "DS2", "DS3"], "bonds" => [0.4025, 0.4765, 0.4955, 0.6335, 0.7005], "phosp" => 2, "apl" => 0.4884),
  "dspe" => Dict("beads" => ["ETH", "PHO", "MID", "DS1", "DS2", "DS3"], "bonds" => [0.3555, 0.4805, 0.4995, 0.6345, 0.7005], "phosp" => 2, "apl" => 0.4368),
  "dspg" => Dict("beads" => ["GLY", "PHO", "MID", "DS1", "DS2", "DS3"], "bonds" => [0.3805, 0.4735, 0.4975, 0.6345, 0.7005], "phosp" => 2, "apl" => 0.4910),
  "dsps" => Dict("beads" => ["SER", "PHO", "MID", "DS1", "DS2", "DS3"], "bonds" => [0.4435, 0.4755, 0.4975, 0.6345, 0.7005], "phosp" => 2, "apl" => 0.4653),

  "popa" => Dict("beads" => [       "PHA", "MID", "PO1", "PO2", "PO3"], "bonds" => [        0.4665, 0.4685, 0.5415, 0.5185], "phosp" => 1, "apl" => 0.5788),
  "popc" => Dict("beads" => ["CHO", "PHO", "MID", "PO1", "PO2", "PO3"], "bonds" => [0.4025, 0.4705, 0.4645, 0.5355, 0.5015], "phosp" => 2, "apl" => 0.6420),
  "pope" => Dict("beads" => ["ETH", "PHO", "MID", "PO1", "PO2", "PO3"], "bonds" => [0.3565, 0.4715, 0.4725, 0.5475, 0.5245], "phosp" => 2, "apl" => 0.5611),
  "popg" => Dict("beads" => ["GLY", "PHO", "MID", "PO1", "PO2", "PO3"], "bonds" => [0.3815, 0.4705, 0.4635, 0.5265, 0.4945], "phosp" => 2, "apl" => 0.6726),
  "pops" => Dict("beads" => ["SER", "PHO", "MID", "PO1", "PO2", "PO3"], "bonds" => [0.4425, 0.4665, 0.4695, 0.5415, 0.5185], "phosp" => 2, "apl" => 0.5770),

  "sopa" => Dict("beads" => [       "PHA", "MID", "SO1", "SO2", "SO3"], "bonds" => [        0.4675, 0.4695, 0.5455, 0.5885], "phosp" => 1, "apl" => 0.5714),
  "sopc" => Dict("beads" => ["CHO", "PHO", "MID", "SO1", "SO2", "SO3"], "bonds" => [0.4025, 0.4715, 0.4645, 0.5365, 0.5745], "phosp" => 2, "apl" => 0.6396),
  "sope" => Dict("beads" => ["ETH", "PHO", "MID", "SO1", "SO2", "SO3"], "bonds" => [0.3565, 0.4725, 0.4705, 0.5485, 0.5915], "phosp" => 2, "apl" => 0.5541),
  "sopg" => Dict("beads" => ["GLY", "PHO", "MID", "SO1", "SO2", "SO3"], "bonds" => [0.3815, 0.4705, 0.4635, 0.5295, 0.5635], "phosp" => 2, "apl" => 0.6684),
  "sops" => Dict("beads" => ["SER", "PHO", "MID", "SO1", "SO2", "SO3"], "bonds" => [0.4435, 0.4675, 0.4715, 0.5445, 0.5845], "phosp" => 2, "apl" => 0.5747),
)

# Argument parser
function parse_commandline()
  s = ArgParseSettings()

  @add_arg_table s begin
    "input"
    arg_type = String
    help = "Input file in TOML format."
    metavar = "INPUT"
    required = true
  end

  return parse_args(s)
end

# Composition parser
function parse_composition(s::String)
  if occursin("/", s)
    a, b = split(s, "/")
    a = parse(Int64, a)
    b = parse(Int64, b)
    n = a // b
  elseif occursin("%", s)
    s = replace(s, "%" => "")
    n = parse(Float64, s) / 100.0
    n = rationalize(n)
  else
    n = parse(Float64, s)
    n = rationalize(n)
  end

  return n
end

# Normalize compositions
function normalize_compositions!(dict::Dict{String,Any})
  total = sum(values(dict))

  for key in keys(dict)
      dict[key] = dict[key] / total
  end
end

# Area per lipid setup
function update_apl(input::Dict{String,Any})
  # Replace values with the ones in the input file
  if haskey(input, "properties")
    for (key, _) in ISOLF
      if haskey(input["properties"], key)
        ISOLF[key]["apl"] = input["properties"][key]["apl"]
      end
    end
  end
end

# Calculate number of lipids
function calc_nlipids(composition::Dict{String,Any}, lx::Float64, ly::Float64)
  # Calculate mean area per lipid
  apl = 0.0
  for (key, val) in composition
    apl += val * ISOLF[key]["apl"]
  end

  # Calculate mean side
  l = sqrt(apl)

  # Calculate number of lipids rounded down
  nx = floor(Int64, lx / l)
  ny = floor(Int64, ly / l)

  # Calculte number of lipids
  n = nx * ny

  return n, nx, ny
end

# Calculate number of lipids
function calc_area(composition::Dict{String,Any}, nlipids::Int64)
  # Calculate total area
  area = 0.0
  for (key, val) in composition
    area += val * ISOLF[key]["apl"] * nlipids
  end
  area = convert(Float64, area)

  return area
end

# Calculate partition of nlipids based on porportions
function calc_partition(composition::Dict{String,Any}, nlipids::Int64)
  # Build buckets
  l = length(composition)
  z = zeros(Rational, l)
  k = keys(composition)
  b = Dict(zip(k, z))
  p = Dict(zip(k, z))

  # Fill buckets
  for i in 1:nlipids
    # Find unhappiest bucket
    kmin = ""
    vmin = 1
    for (key, val) in composition
      p[key] -= val
      if p[key] < vmin
        vmin = p[key]
        kmin = key
      end
    end
    # Fill unhappiest bucket
    b[kmin] += 1
    # Update proportions
    for key in keys(composition)
      p[key] = b[key] // i
    end
  end

  # Build result
  partition = Dict()
  for (key, val) in b
    partition[key] = convert(Int64, val)
  end

  return partition
end

# Calculate positions
function calc_particles_positions(composition::Dict{String,Any}, nlipids::Int64, nx::Int64, ny::Int64, lx::Float64, ly::Float64, aid::Int64, mid::Int64, invert::Bool)
  # Calculate partition
  partition = calc_partition(composition, nlipids)
  # Create array of lipids
  lipids = String[]
  for (key, val) in partition
    lipids = vcat(lipids, repeat([key], val))
  end
  # Randomize array
  lipids = shuffle(lipids)
  # Calculate number of particles
  np = 0
  for (key, val) in partition
    np += length(ISOLF[key]["beads"]) * val
  end
  # Calculate sides
  sx = lx / nx
  sy = ly / ny
  # Calculate particle positions
  p = Array{Particle,1}(undef, np)
  m = 0
  for i in 0:(nx-1)
    for j in 0:(ny-1)
      # Update lipid number
      mid += 1
      # Calculate lipid id
      lid = i * ny + j + 1
      key = lipids[lid]
      # Calculate xy position
      x = i * sx
      y = j * sy
      # Get number of beads
      n = length(ISOLF[key]["beads"])
      # Create particles
      lname = key
      z = 0.0
      for k in 1:(ISOLF[key]["phosp"]-1)
        z += ISOLF[key]["bonds"][k]
      end
      for k in 1:n
        # Update particle number
        aid += 1
        # Build particle
        m += 1
        bname = ISOLF[key]["beads"][k]
        p[m] = Particle(lname, bname, mid, aid, x, y, z)
        if k < n
          z -= ISOLF[key]["bonds"][k]
        end
      end
    end
  end
  # Shift positions to center in the origin
  zmin = 0.0
  for particle in p
    z = particle.z
    zmin = z < zmin ? z : zmin
  end
  gap = min(sx, sy) / 2.0
  xsh = (sx - lx) / 2.0
  ysh = (sy - ly) / 2.0
  zsh = gap - zmin
  fac = invert ? -1.0 : 1.0
  for i in 1:np
    p[i].x += xsh
    p[i].y += ysh
    p[i].z += zsh
    p[i].z *= fac
  end

  return p, lipids, aid, mid
end

# GRO writer
function write_gro(file_name::String, system_name::String, particles::Vector{Particle}, l::Vector{String})
  # Calculate box size
  xmin = 0.0
  xmax = 0.0
  ymin = 0.0
  ymax = 0.0
  zmin = 0.0
  zmax = 0.0
  for p in particles
    xmin = p.x < xmin ? p.x : xmin
    xmax = p.x > xmax ? p.x : xmax

    ymin = p.y < ymin ? p.y : ymin
    ymax = p.y > ymax ? p.y : ymax

    zmin = p.z < zmin ? p.z : zmin
    zmax = p.z > zmax ? p.z : zmax
  end
  bx = xmax - xmin
  by = ymax - ymin
  bz = zmax - zmin
  # Get lipid names
  keys = unique(l)
  # Write gro file
  new_mid = 0
  act_mid = 0
  act_aid = 0
  open(file_name, "w") do io
    @printf(io, "%s, t = %16.3f \n", system_name, 0)
    @printf(io, "%12d \n", length(particles))
    for key in keys
      for p in particles
        if p.mol != key
          continue
        end
        if new_mid != p.mid
          new_mid = p.mid
          act_mid = act_mid + 1
        end
        act_aid = act_aid + 1
        @printf(io, "%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f \n",
          act_mid % 100000,
          p.mol,
          p.atm,
          act_aid % 100000,
          p.x,
          p.y,
          p.z,
          0.0,
          0.0,
          0.0
        )
      end
    end
    @printf(io, "%15.4f%15.4f%15.4f \n\n", bx, by, bz)
  end
end

# PDB writer
function write_pdb(file_name::String, system_name::String, particles::Vector{Particle})
  open(file_name, "w") do io
    @printf(io, "TITLE     %70s\n", rpad(system_name, 70, " "))
    for p in particles
      @printf(io,
        "ATOM  %5d %4s%1s%4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%10s%2s%2s \n",
        p.aid % 100000,
        p.atm,
        " ",
        rpad(uppercase(p.mol), 4, " "),
        "A",
        p.mid % 10000,
        " ",
        p.x * 10,
        p.y * 10,
        p.z * 10,
        0.0,
        0.0,
        "",
        "",
        "")
    end
    @printf(io, "%s", "TER\n")
  end
end

# TOP writer
function write_top(file_name::String, system_name::String, l::Vector{String})
  # Count lipid numbers
  keys = unique(l)
  lips = Dict{String,Int64}()
  for key in keys
    lips[key] = 0
  end
  for s in l
    lips[s] = lips[s] + 1
  end

  # Write top file
  open(file_name, "w") do io
    print(io, "; Implicit Solvent Lipid Forcefield (iSoLF)\n")
    print(io, "#include \"./param/isolf.itp\"\n")
    print(io, "\n")

    print(io, "[ system ]\n")
    print(io, system_name * "\n")
    print(io, "\n")

    print(io, "[ molecules ]\n")
    for (key, val) in lips
        @printf(io, "%s %d\n", uppercase(key), val)
    end
    print(io, "\n")

    print(io, "[ cg_ele_chain_pairs ]\n")
    @printf(io, "ON 1 - %d : 1 - %d\n", length(l), length(l))
  end
end

# iSoLF membrane builder
function build_membrane(input::Dict{String,Any})
  # Parse upper layer composition into rationals
  upper = input["membrane"]["upper"]["composition"]
  for (key, val) in upper
    upper[key] = parse_composition(val)
  end
  normalize_compositions!(upper)

  # Parse lower layer composition into rationals
  lower = input["membrane"]["lower"]["composition"]
  for (key, val) in lower
    lower[key] = parse_composition(val)
  end
  normalize_compositions!(lower)

  # Parse membrane lengths
  if haskey(input["membrane"], "size")
    lx = get(input["membrane"]["size"], "x", 0.0)
    ly = get(input["membrane"]["size"], "y", 0.0)
  else
    lx = 0.0
    ly = 0.0
  end

  # Parse number of lipids
  nu = get(input["membrane"]["upper"], "nlipids_per_side", 0)
  nl = get(input["membrane"]["lower"], "nlipids_per_side", 0)
  nu = nu^2
  nl = nl^2

  # Parse output
  out_path = get(input["membrane"]["output"], "path", "./")
  out_name = get(input["membrane"]["output"], "name", "out")
  out_fmts = get(input["membrane"]["output"], "formats", ["gro", "top"])

  # Process paths
  mkpath(out_path)
  out_files = Dict{String,String}()
  for fmt in out_fmts
    out_files[fmt] = joinpath(out_path, out_name * "." * fmt)
  end

  # Parse are per lipid
  update_apl(input)

  # Calculate membrane lengths and number of lipids if necessary
  if lx * ly > 0.0 && nu * nl == 0
    nu, nux, nuy = calc_nlipids(upper, lx, ly)
    nl, nlx, nly = calc_nlipids(lower, lx, ly)
  elseif lx * ly == 0.0 && nu * nl > 0
    au = calc_area(upper, nu)
    al = calc_area(lower, nl)
    area = max(au, al)
    lx = sqrt(area)
    ly = sqrt(area)
    lu = sqrt(au / nu)
    ll = sqrt(al / nl)
    nux = floor(Int64, lx / lu)
    nuy = floor(Int64, ly / lu)
    nlx = floor(Int64, lx / ll)
    nly = floor(Int64, ly / ll)
  end

  # Calculte positions for each lipid layer
  aid = 0
  mid = 0
  pu, lipu, aid, mid = calc_particles_positions(upper, nu, nux, nuy, lx, ly, aid, mid, false)
  pl, lipl, aid, mid = calc_particles_positions(lower, nl, nlx, nly, lx, ly, aid, mid, true)

  # Concatenat layers
  p = vcat(pu, pl)
  l = vcat(lipu, lipl)

  # Write files
  for (key, val) in out_files
    if key == "gro"
      write_gro(val, "CG membrane model", p, l)
    end
    if key == "pdb"
      write_pdb(val, "CG membrane model", p)
    end
    if key == "top"
      write_top(val, "CG membrane model", l)
    end
  end
end

function main()
  # Parse arguments
  args = parse_commandline()

  # Open input file
  io = open(args["input"], "r")

  # Load content into a string
  content = read(io, String)

  # Parse content as a TOML file
  input = TOML.parse(content)

  # Build membrane is required
  if haskey(input, "membrane")
    build_membrane(input)
  end
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
