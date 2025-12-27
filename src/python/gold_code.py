import numpy as np

class LFSR:
    def __init__(self, poly, initial_state):
        """
        Fibonacci LFSR.
        poly: List of exponents for the polynomial terms.
              e.g. [10, 3, 0] -> x^10 + x^3 + 1
              The '0' implies the input +1, and the highest degree is the length.
        initial_state: List/Array of 0s and 1s, length = degree.
        """
        self.degree = poly[0]
        # Taps: The powers that appear in the polynomial (excluding 0 and N if implicit logic used, 
        # but usually we need the intermediate taps).
        # For [10, 3, 0], we tap the output (which is always used in feedback) and the bit corresponding to x^3.
        # In a shift register [0..N-1], if we shift right, the output is index 0.
        # If we shift left, output is index N-1.
        # Let's assume standard shift-left for bitwise convenience (or right).
        
        # We'll use an integer mask implementation for speed.
        self.taps = [p for p in poly if 0 < p < self.degree]
        
        # Convert state to integer
        # init_state = [b0, b1, ..., bN-1]
        self.state = 0
        for b in initial_state:
            self.state = (self.state << 1) | int(b)
        
        self.mask = 0
        for t in self.taps:
            # If poly is x^10 + x^3 + 1. 
            # In a 10-bit register usually corresponding to x^1...x^10?
            # Feedback = bit_at_pos(10) XOR bit_at_pos(3).
            # But we only store 10 bits.
            # Usually: Feedback = State[Degree-Tap]. 
            # Let's simpler logic: 
            # "One-to-many" or "Many-to-one"? Fibonacci is Many-to-one.
            # Feedback bit = XOR sum of tapped stages.
            # Stages are usually 1..N.
            # Tap 3 means stage 3.
            # If we store state as integer, bit 0 is stage 1?
            self.mask |= (1 << (t - 1))     

    def step(self):
        # Read the output bit (usually the MSB or LSB depending on shift direction)
        # Let's assume shift LEFT. MSB (bit degree-1) falls off.
        out_bit = (self.state >> (self.degree - 1)) & 1
        
        # Tapped bits:
        # We XOR the bits at tap positions.
        # But wait, usually the "output" is fed back xor'd with other taps.
        # Feedback = out_bit XOR (other taps).
        
        # Actually, simpler:
        # Feedback = XOR( bits at poly-taps ).
        # For x^10 + x^3 + 1, we tap stage 10 (output) and stage 3.
        # feedback = (state >> (10-1)) ^ (state >> (3-1))
        
        feedback = out_bit
        # XOR with other taps
        for t in self.taps:
            tap_bit = (self.state >> (t - 1)) & 1
            feedback ^= tap_bit
            
        # Shift and insert feedback at LSB
        self.state = ((self.state << 1) & ((1 << self.degree) - 1)) | feedback
        return out_bit

    def generate(self, n):
        out = np.zeros(n, dtype=np.uint8)
        for i in range(n):
            out[i] = self.step()
        return out

class GoldCodeGenerator:
    def __init__(self, poly1, state1, poly2, state2):
        self.lfsr1 = LFSR(poly1, state1)
        self.lfsr2 = LFSR(poly2, state2)
        
    def generate(self, length):
        seq1 = self.lfsr1.generate(length)
        seq2 = self.lfsr2.generate(length)
        return np.bitwise_xor(seq1, seq2)
