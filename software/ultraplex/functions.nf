def initOptions(Map args) {
    def Map options = [:]
    options.args            = args.args ?: ''
    options.publish_dir     = args.publish_dir ?: 'ultraplex'
    options.publish_results = args.publish_results ?: 'all'
    return options
}
