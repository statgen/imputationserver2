import java.security.SecureRandom

class PasswordCreator {

    static SecureRandom rng = new SecureRandom()

    static String createPassword() {
        return createPassword(6, 6, 3, 3)
    }

    static String createPassword(int uppercaseLetters, int lowercaseLetters, int numbers, int symbols) {
        List<String> password = []

        password += take('ABCDEFGHIJKLMNOPQRSTUVWXYZ', uppercaseLetters)
        password += take('abcdefghijklmnopqrstuvwxyz', lowercaseLetters)
        password += take('0123456789', numbers)
        password += take('!@#%&*()_+[]{};:.<>?-=', symbols)

        Collections.shuffle(password, rng)
        return password.join()
    }

    private static List<String> take(String population, int n) {
        List<String> options = population.toList()
        Collections.shuffle(options, rng)
        return options.take(n)
    }

}
