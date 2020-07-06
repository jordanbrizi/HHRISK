// This file is required by the index.html file and will
// be executed in the renderer process for that window.
// All of the Node.js APIs are available in this process.

const pkg = () => require('./package.json')
const tabela = () => require('./tabela.json')
const defaultLang = Intl.DateTimeFormat().resolvedOptions().locale
const lang = (attr = defaultLang) => require(`./lang/${attr}.json`)

module.exports = {pkg, tabela, lang}