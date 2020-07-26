const { ipcRenderer, remote, shell } = require('electron')
const defaultLang = Intl.DateTimeFormat().resolvedOptions().locale
const pkg = () => require('../package')
const resultsPath = arg => remote.app.getAppPath() + '\\bin\\Results\\' + arg
const mainPath = arg => remote.app.getAppPath() + '\\' + arg

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

const lang = (attr = defaultLang) => {
	if(attr != null) {
		return require(mainPath(`ui\\lang\\${attr}`))
	}
	else {
		return require(mainPath(`ui\\lang\\${defaultLang}`))	
	}
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

module.exports = {pkg, lang, ipcRenderer, shell, remote}