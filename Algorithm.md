```plaintext
Algorithm 1  The main framework of the HEA-ELITE algorithm
    Input: A Latin square graph G = (V, E), Gene pool size M, Mutation parameter K
    Output: A coloring solution C

    The best found solution C* <- Initialize(G)
    Gene pool S <- {S_i <- Initialize(G)}
    Elite pool E <- EmptySet

    while elapsed time < time limit do
        Select the best sample X <- min(S)
        Select a random sample Y from S
        C <- Genetic(X, Y) // Algorithm 2
        C <- TabuSearch(C)

        if |CE(C)| = 0 then
            return C

        if |CE(C)| < |CE(C*)| then
            C* <- C

        S <- S \cup {C}
        if |S| = 2M then
            Keep the best M samples of S
            if min(S) = max(S) then
                Select a random sample e from S
                E <- E \cup {e}
                if |E| = M:
                    S = E
                    clear E

        Randomly select K samples for mutation


Algorithm 2  Genetic algorithm
    Input: coloring solutions X and Y
    Output: A coloring solution Z

    Conflict position set S <- EmptySet
    for i = 1 to n do
        Z_i = X_i
        if \delta|CE(Z)| > 0 then
            S <- S \cup {i}
    
    for all i \in S do
        Z_i = Y_i
    return Z
```