# =====================================================================
# === Structures
  
  struct ProblemConstants{T}
    MJDᵢ::T
    g₀::T
    μ::T
    R::T
    Io::PlanetaryConstantsModel
    Europa::PlanetaryConstantsModel
    Ganymede::PlanetaryConstantsModel
    Callisto::PlanetaryConstantsModel
  end
  
  # =====================================================================
  # === Constructors
  
  function import_constants()
    # Values from GTOC6 overview
    μⱼ = 126686534.92180
    tᵢ = 58849.0
    g = 9.80665
    R = 71492.0
  
    # Creating Moons
    a = 422029.68714001
    Io = PlanetaryConstantsModel(
      5959.916, 
      1826.5, 
      sqrt(μⱼ/a^3),
      Keplerian(
        a, 
        4.308524661773E-03, 
        40.11548686966E-03,
        -79.640061742992, 
        37.991267683987,    
        286.85240405645, 
        0.0)
      )
  
    a = 671224.23712681
    Europa = PlanetaryConstantsModel(
      3202.739, 
      1561.0, 
      sqrt(μⱼ/a^3),
      Keplerian(
        a,   
        9.384699662601E-03, 
        0.46530284284480,
        -132.15817268686, 
        -79.571640035051, 
        318.00776678240, 
        0.0)
      )
  
    a = 1070587.4692374
    Ganymede = PlanetaryConstantsModel(
      9887.834, 
      2634.0, 
      sqrt(μⱼ/a^3),
      Keplerian(
        a, 
        1.953365822716E-03, 
        0.13543966756582,
        -50.793372416917, 
        -42.876495018307, 
        220.59841030407, 
        0.0)
      )
  
    a = 1883136.6167305
    Callisto = PlanetaryConstantsModel(
      7179.289, 
      2408.0, 
      sqrt(μⱼ/a^3),
      Keplerian(
        a, 
        7.337063799028E-03, 
        0.25354332731555,
        86.723916616548, 
        -160.76003434076, 
        321.07650614246, 
        0.0)
      )
  
  
    return ProblemConstants(tᵢ, g, μⱼ, R, 
      Io, Europa, Ganymede, Callisto)
  end