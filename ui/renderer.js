const { ipcRenderer } = require('electron')
const path = require('path')
const pkg = () => require('./package')
const Exposure = () => require(path.resolve(`../bin/Results/Exposure`))
const defaultLang = Intl.DateTimeFormat().resolvedOptions().locale

const lang = (attr = defaultLang) => {
	if(attr != null) {
		return require(`./lang/${attr}`)
	}
	else {
		return require(`./lang/${defaultLang}`)	
	}
}
module.exports = {pkg, lang, Exposure, ipcRenderer}