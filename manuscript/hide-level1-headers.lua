-- hide-level1-headers.lua
-- Remove all level-1 headings from the output.

function Header(el)
  if el.level == 1 then
    -- returning {} deletes this header node
    return {}
  end
  return el
end
