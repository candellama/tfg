import oqs
import time

msg = b"prueba de firma digital con dilithium libreria oqs tfg"

# vector of security levels to test
security_levels = ["ML-DSA-44", "ML-DSA-65", "ML-DSA-87"]
# initialize results dictionary
results = {}

# loop testing each security level
for level in security_levels:

    # generate key pair and signature of msg
    with oqs.Signature(level) as signer:
        pk = signer.generate_keypair()
        time_start = time.time()
        sig = signer.sign(msg)
        signer.last_sign_time = time.time() - time_start

    # verify signature
    with oqs.Signature(level) as verifier:
        ok = verifier.verify(msg, sig, pk)

    results[level] = {
        "signature_valid": ok,
        "key_size": len(pk),
        "signature_size": len(sig),
        "signing_time": signer.last_sign_time
    }

print("Results for ML-DSA (Dilithium) signature scheme:")
for level, data in results.items():
    print(f"Security Level: {level}")
    print(f"  Valid signature: {data['signature_valid']}")
    print(f"  Public Key Size: {data['key_size']} bytes")
    print(f"  Signature Size: {data['signature_size']} bytes")
    print(f"  Time to Sign: {data['signing_time']} seconds")

