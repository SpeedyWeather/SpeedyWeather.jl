// .vitepress/theme/index.ts
import { h } from 'vue'
import DefaultTheme from 'vitepress/theme'
import type { Theme as ThemeConfig } from 'vitepress'
import 'virtual:mathjax-styles.css';

import { 
  NolebaseEnhancedReadabilitiesMenu, 
  NolebaseEnhancedReadabilitiesScreenMenu, 
} from '@nolebase/vitepress-plugin-enhanced-readabilities/client'

import { enhanceAppWithTabs } from 'vitepress-plugin-tabs/client'

import '@nolebase/vitepress-plugin-enhanced-readabilities/client/style.css'
import './style.css'
import './docstrings.css'

export const Theme: ThemeConfig = {
  extends: DefaultTheme,
  Layout() {
    return h(DefaultTheme.Layout, null, {
      'nav-bar-content-after': () => [
        h(NolebaseEnhancedReadabilitiesMenu),
      ],
      'nav-screen-content-after': () => h(NolebaseEnhancedReadabilitiesScreenMenu),
    })
  },
  enhanceApp({ app, router, siteData }) {
    enhanceAppWithTabs(app);
  }
}
export default Theme
