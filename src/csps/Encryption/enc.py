import Crypto
from Crypto.PublicKey import RSA
import ast

f = open('public.pem', 'rb')
publickey = RSA.importKey(f.read())
f.close()

encrypted = publickey.encrypt('AKIAI7XOEFVYPOUQKGYQ', 32)
print (encrypted)

f2 = open ('aws_key.txt', 'w')
f2.write(str(encrypted))
f2.close()

# Testing decryption

passphrase = str(raw_input("Please enter the passphrase for the private key: "))

f1 = open('private.pem', 'rb')
privatekey = RSA.importKey(f1.read(), passphrase=passphrase)
f1.close()

f3 = open('aws_key.txt', 'r')
message = f3.read()
f3.close()

decrypted = privatekey.decrypt(ast.literal_eval(str(message)))
print (decrypted)
