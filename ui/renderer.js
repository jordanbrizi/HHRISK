const { ipcRenderer, remote, shell } = require('electron')
const defaultLang = Intl.DateTimeFormat().resolvedOptions().locale

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

module.exports = { defaultLang, ipcRenderer, shell, remote }