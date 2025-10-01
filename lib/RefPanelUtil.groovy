import org.yaml.snakeyaml.Yaml

class RefPanelUtil {

    public static Map<String, Object> loadFromFile(filename) {
        println "Loading reference panel from file ${filename}..."
        def params_file = new File(filename)
        def parameter_yaml = new FileInputStream(params_file)
        def panel = new Yaml().load(parameter_yaml)

        HashMap<String, String> environment = [:]
        def folder = params_file.getParentFile().getAbsolutePath()
        environment.put('CLOUDGENE_APP_LOCATION', folder)

        RefPanelUtil.resolveEnv(panel.properties, environment)
        return panel.properties
    }

    public static void resolveEnv(properties, environment) {
        for (String property : properties.keySet()) {
            Object propertyValue = properties.get(property)
            if (propertyValue instanceof String) {
                propertyValue = RefPanelUtil.env(propertyValue.toString(), environment)
            } else if (propertyValue instanceof Map) {
                resolveEnv(propertyValue, environment)
            }
            properties.put(property, propertyValue)
        }
    }

    public static String env(String value, Map<String, String> variables) {
        String updatedValue = value + ''
        for (String key : variables.keySet()) {
            updatedValue = updatedValue.replaceAll('\\$\\{' + key + '\\}', variables.get(key))
        }
        return updatedValue
    }

}
