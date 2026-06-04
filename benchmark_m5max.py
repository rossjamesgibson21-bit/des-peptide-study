"""
E1.1 — M5 Max throughput benchmark (v2)
Loads from inpcrd (with box vectors), minimises, then benchmarks.
"""
import time
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import os

SYSTEMS = {
    "GGE_water":      ("systems/GGE_water/GGE_water.prmtop",        "systems/GGE_water/GGE_water.inpcrd"),
    "CME_reline":     ("systems/CME_reline/CME_reline.prmtop",      "systems/CME_reline/CME_reline.inpcrd"),
    "YIY_glyceline":  ("systems/YIY_glyceline/YIY_glyceline.prmtop", "systems/YIY_glyceline/YIY_glyceline.inpcrd"),
}

PLATFORMS = ["CPU", "OpenCL"]
STEPS = 5000
TIMESTEP = 2.0
WARMUP_STEPS = 500

os.chdir(os.path.expanduser("~/des-peptide-study"))

results = []

for sys_name, (prmtop_path, inpcrd_path) in SYSTEMS.items():
    print(f"\n{'='*60}")
    print(f"System: {sys_name}")
    
    prmtop = app.AmberPrmtopFile(prmtop_path)
    inpcrd = app.AmberInpcrdFile(inpcrd_path)
    n_atoms = prmtop.topology.getNumAtoms()
    print(f"  Atoms: {n_atoms}")
    print(f"  Box vectors: {inpcrd.boxVectors is not None}")
    
    for plat_name in PLATFORMS:
        print(f"  Platform: {plat_name} ... ", end="", flush=True)
        
        try:
            platform = mm.Platform.getPlatformByName(plat_name)
            
            system = prmtop.createSystem(
                nonbondedMethod=app.PME,
                nonbondedCutoff=1.0 * unit.nanometers,
                constraints=app.HBonds,
            )
            system.addForce(mm.MonteCarloBarostat(1.0 * unit.bar, 300 * unit.kelvin))
            
            integrator = mm.LangevinIntegrator(
                300 * unit.kelvin,
                1.0 / unit.picoseconds,
                TIMESTEP * unit.femtoseconds,
            )
            
            simulation = app.Simulation(prmtop.topology, system, integrator, platform)
            simulation.context.setPositions(inpcrd.positions)
            if inpcrd.boxVectors is not None:
                simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
            
            # Minimise to remove any clashes
            print("minimising... ", end="", flush=True)
            simulation.minimizeEnergy(maxIterations=1000)
            
            # Warmup (not timed)
            simulation.step(WARMUP_STEPS)
            
            # Timed run
            start = time.time()
            simulation.step(STEPS)
            elapsed = time.time() - start
            
            ns_simulated = STEPS * TIMESTEP * 1e-6
            ns_per_day = ns_simulated / elapsed * 86400
            
            print(f"{ns_per_day:.1f} ns/day  ({elapsed:.1f} s for {ns_simulated:.3f} ns)")
            results.append((sys_name, n_atoms, plat_name, ns_per_day, elapsed))
            
            del simulation, integrator, system
            
        except Exception as e:
            print(f"FAILED: {e}")
            results.append((sys_name, n_atoms, plat_name, 0, 0))

print(f"\n{'='*60}")
print(f"{'System':<20} {'Atoms':>6} {'Platform':<10} {'ns/day':>10}")
print(f"{'-'*60}")
for sys_name, n_atoms, plat_name, ns_day, elapsed in results:
    print(f"{sys_name:<20} {n_atoms:>6} {plat_name:<10} {ns_day:>10.1f}")
print(f"{'='*60}")
