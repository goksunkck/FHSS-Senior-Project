# Future Attack Strategy: Neural Peeling (Hypothesis)

## Concept
The user proposed a "Divide and Conquer" strategy:
1. Use the Neural Network to learn partial information (Parity).
2. Use this information to mathematically isolate ("Peel") one of the constituent LFSRs.
3. Train a new model on the remainder.

## Mathematical Basis
A Gold Code $G(t)$ is formed by:
$$ G(t) = L_1(t) \oplus L_2(t) $$

Our Experiment 6 showed the Transformer reached 25% accuracy by learning the **Parity Behavior** of $G(t)$. 
If this Parity Prediction $P_{pred}(t)$ has a correlation with $L_1(t)$:
$$ Corr(P_{pred}, L_1) = \epsilon > 0 $$

Then we can treat $L_2(t)$ as "noise" and solve for $L_1(t)$ using a Fast Correlation Attack.

## The Algorithm
1. **Inference:** Run the trained Transformer on the Gold Code signal. Record the "Soft Predictions" (logits) before the final decision.
2. **Correlation:** Cross-correlate these soft predictions with all possible phases of the known generator polynomial for LFSR 1 ($Poly_A$).
   - Complexity: $2^{10}$ (Linear search for Phase).
3. **Identification:** The correct phase of $L_1$ will verify with a correlation spike.
4. **Peeling (XOR Subtraction):**
   $$ \text{Residual} = \text{Signal} \oplus \text{Reconstructed } L_1 $$
   $$ \text{Residual} = (L_1 \oplus L_2) \oplus L_1 = L_2 $$
5. **Final Solve:** The Residual is now a pure M-Sequence ($L_2$).
   - Feed this Residual into the Transformer.
   - Based on Experiment 5, the model will solve this with **99.9% accuracy**.

## User Question: Direct Extraction?
**Question:** "Why try every correlation? Can't we just extract the LFSR the network found and subtract it?"

**Answer (The Symmetry Problem):** 
The difficulty lies in the fact that the Neural Network learns the **Sum** ($L_1 \oplus L_2$), not the individual components.
- Standard Gold Codes use two LFSRs of equal length and strength.
- The Network cannot distinguish whether a "1" came from $L_1$ or $L_2$.
- Therefore, the Network's prediction is: $Pred \approx L_1 \oplus L_2$.
- Attempting to subtract this prediction directly:
  $$ Input - Pred \approx (L_1 \oplus L_2) - (L_1 \oplus L_2) \approx 0 $$
  This simply cancels the signal.

**Exception (Asymmetric Codes):**
The user's proposed strategy **would** work if the system was asymmetric (e.g. $L_1$ is much shorter/simpler than $L_2$).
- The Network would "latch on" to the simpler pattern ($L_1$) and treat $L_2$ as noise.
- Then $Pred \approx L_1$.
- Peeling: $(L_1 \oplus L_2) - Pred \approx L_2$.
- This effectively "strips" the weak code to reveal the strong code.

## Conclusion
The "Partial Break" (25%) is not a dead end. It is a sufficient statistic to enable a full cryptographic break using a Neural-Algebraic hybrid approach, reducing the search space from $2^{20}$ (Unbreakable) to $2^{10} + 2^{10}$ (Trivial).
