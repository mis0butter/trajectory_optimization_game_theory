# trajectory_optimization_game_theory
For game theory class 

### Schedule 

Bullet points are goals to complete by the date 

Oct 23: 
- simulator in Julia with propagator 
- 2-body (J2 and J3) 
- model chaser + evader 
- add fuel/maneuver limits 
- min distance between satellites 

... 
- experiment with plotting tools at some point 


Nov 27 
- code done 

Nov 28 
- presentation done 

Nov 29
- present 

## Oct 24, 2023 

Static vs. dynamic game ? 
- Static game uses a-e-i and some q-Q-something thing to compute reachability sets for evader / pursuer 
  - Reachability sets help determine if evasion or capture *can* be guaranteed (but still up to user to compute optimal control sequence) 
- Dynamic game 
  - We just propagate dynamics 
  - Still need to formulate cost function: mass-optimal? time-optimal? 

Bryn technical advice: 
- This is how you formulate cost function for single player optimal control:   
$$ J = \sum_{i=1}^{N} \Delta v_i^2 - K_{terminal} $$ 
- Constraints will be dynamics  

Path forward: 
- Static game where we solve above cost function for delta v at each segment  
- $K_{terminal}$ is basically $K_{distance}$ 
- Evader (P1) wants to minimize $J$ --> $ \Delta v $ and maximize $ K_{distance} $, thus minus sign 
- Pursuer (P2) wants to maximize $J$ 
- terminal constraint: 
$$ \Delta \bar{r} = \bar{0} = \bar{r}_E - \bar{r}_P $$ 

Questions: 
- Constraint dynamics will have to be integrated for each time segment and then $\Delta v$ applied between each segment 
$$ c_1 = \int_{t_0}^{t_1} f( \Delta v_1 )dt + ... + \int_{t_f-1}^{t_f} f( \Delta v_f ) dt $$ 
  - How to implement this kind of constraint? 
- Open loop? Feedback? 
- Perfect information at each state? 



## Oct 17, 2023 

### Actions 

Sofia: 
- simulator in Julia with propagator 
- 2-body (J2 and J3) 
- model chaser + evader 

Junette: 
- add fuel/maneuver limits 
- min distance between satellites 

### Notes 

The point of the class is: 
- non-cooperative game theory 
- minimize different objective functions 
  - based on fuel 

Minimum success criteria for simulator dynamics:
- IC target 
- IC chaser 
- propagator 
- 2-body dynamics 
- fuel/maneuver limits 
- min distance between sets 
- shadow (observability) 

Nice to have: 
- shadow model (observability) 
- atmosphere (cannonball) 








