#[cfg(test)]
mod tests {
    use crate::preprocessed_inputs::preprocessed_inputs;
    use crate::prover::prover;
    use crate::verifier::verifier;

    #[test]
    fn test() {
        let preprocessed_inputs = preprocessed_inputs();
        let proof = prover(&preprocessed_inputs);
        let result = verifier(preprocessed_inputs, proof);

        if result {
            println!("{}", "The proof is accepted!");
        } else {
            println!("{}", "The proof is not accepted!");
        }
    }
}
